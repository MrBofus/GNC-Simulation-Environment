#```````````````````````````````````````````````````````````````````````````````````````````````````````````````##
# library imports
#
# --> removed pyatmos, doesn't work with python 3.9.13 (?)

from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Earth
from astropy import units as u
# from pyatmos import download_sw_jb2008, read_sw_jb2008
# from pyatmos import jb2008
import numpy as np
import pyproj
import random

# swfile = download_sw_jb2008() 

#```````````````````````````````````````````````````````````````````````````````````````````````````````````````##
# helper functions
#

# compute the magnitude of a vector
def magnitude(vec):

    # initialize magnitude to zero
    mag = 0

    # iterate through vector and sum the squares of the components
    for v in vec:
        mag += v**2
    
    # return the magnitude of the vector
    return np.sqrt(mag)


# compute a unit vector along a basis vector
def normalize(vec):
    
    # determine magnitude of basis vector
    m = magnitude(vec)
    
    # initialize unit vector
    v_ = []

    # iterate through basis vector and divide each component
    # by the magnitude
    for v in vec:
        v_.append(v/m)
    
    # return unit vector
    return v_


# transform between ECI and ECEF coordinates
def ECI_to_ECEF(r, t):
    
    # account for rate of Earth's rotation
    gamma = (360/(23.9345*3600)) * t
    gamma *= np.pi/180
    
    # perform the rotation using rotation matrix
    x = r[0]*np.cos(gamma) - r[1]*np.sin(gamma) # meters
    y = r[0]*np.sin(gamma) + r[1]*np.cos(gamma) # meters
    z = r[2]                                    # meters
    
    # define the ellipsoid projection
    transformer = pyproj.Transformer.from_crs({"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
                                              {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'})
    
    # determine ECEF coordinates
    lon, lat, alt = transformer.transform(x, y, z, radians=False)
    
    return lon, lat, alt


# determine the density of air given altitude
def idealGas(height):
   
    # https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html

    # determine temperature and pressure using magic numbers
    T = -131.21 + 0.00299*height
    p = 2.488 * ( (T + 273.1)/216.6 ) ** -11.388

    # determine air density using magic numbers
    rho = p / (0.2869 * (T + 273.1))
    return rho


# compute drag acting on satellite
def compute_drag(r, v, A, m, t):
    
    # read space weather data to determine air density  [depreciated]
    # swdata = read_sw_jb2008(swfile)
    
    # define time and determine ECEF location
    time = '2014-07-22 22:18:45'
    lon, lat, alt = ECI_to_ECEF(r, t)

    # determine air density using space weather data    [depreciated]
    # jb08 = jb2008(time, (lat, lon, alt/10**3), swdata)
    # rho = jb08.rho

    # determine air density using approximation
    rho = idealGas(alt)
    rho += rho * (0.5 - random.random()) / 5

    # drag coefficient
    Cd = 2.2
    
    # return magnitude of drag force acting on satellite
    return (1/2) * rho * Cd * A * v**2


# integrate acceleration
def integrate_accel(accel, velocity, dt):

    # initialize velocity vector
    v_candidate = []

    # iterate through velocity and acceleration
    for i in range(3):
        v_candidate.append(velocity[i] + accel[i]*dt)
    
    # return updated velocity
    return v_candidate



#```````````````````````````````````````````````````````````````````````````````````````````````````````````````##
# force propagator
#

def fpropagate(orbit, f_thrust, mass, A, t, dt):
    orbit = orbit.propagate(dt << u.second)
    
    r, v = (orbit.r << u.meter).value, (orbit.v << u.meter/u.second).value
    
    f_drag_magnitude = compute_drag(r, magnitude(v), A, mass, t)
    f_drag = -1 * f_drag_magnitude * np.array(normalize(v))
    
    f_total = f_drag + f_thrust
    
    v = integrate_accel(f_total/mass, v, dt)
    
    return Orbit.from_vectors(Earth, r << u.meter, v << u.meter/u.second)


if __name__ == "__main__":

    import gnc_library as gnc
    from poliastro.twobody.orbit import Orbit
    from poliastro.bodies import Earth
    from astropy import units as u
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np


    def spiral_in(f_mag, orbit, target):
        altitude = (orbit.a << u.km).value - 6378.1
        if altitude > target:
            thrust = -1 * f_mag * np.array( gnc.normalize(orbit.v.value) )
        
        else:
            thrust = [0, 0, 0]
        
        return thrust

    def spiral_out(f_mag, orbit, target):
        altitude = (orbit.a << u.km).value - 6378.1
        if altitude < target:
            thrust = f_mag * np.array( gnc.normalize(orbit.v.value) )
        
        else:
            thrust = [0, 0, 0]
        
        return thrust


    R_earth = 6378.1

    a_initial = (565.25 + R_earth) << u.km
    e_initial = 10**-4 << u.one
    i_initial = 45 << u.deg
    w_initial = 0 << u.deg
    O_initial = 0 << u.deg
    nu_initial = 0 << u.deg 

    orbit = Orbit.from_classical(Earth, a_initial, e_initial, i_initial, O_initial, w_initial, nu_initial)

    v_initial = gnc.magnitude((orbit.v << u.meter/u.second).value)

    f_mag = 1.3
    mass = 80
    A = 2


    df = pd.DataFrame()

    print('\n')

    xdata, ydata, zdata = [], [], []

    t = 0
    dt = 60*1
    t_max = 8*90*60
    while t < t_max:

        f_thrust = spiral_out(f_mag, orbit, 1500)
        orbit = fpropagate(orbit, f_thrust, mass, A, t, dt)

        tempdf = pd.DataFrame({'time':[t],
                               'x':[(orbit.r[0] << u.km).value],
                               'y':[(orbit.r[1] << u.km).value],
                               'z':[(orbit.r[2] << u.km).value],
                               'altitude':[(orbit.a << u.km).value - R_earth],
                               'inclination':[(orbit.inc << u.deg).value]})
        df = pd.concat([df, tempdf])

        print("\r" + str(int(100*t/t_max)) + "% complete", end='')
        t += dt

    print('\n')
    v_final = gnc.magnitude((orbit.v << u.meter/u.second).value)

    print('v_initial: ' + str(round(v_initial, 1)) + 'm/s')
    print('v_final: ' + str(round(v_final, 1)) + 'm/s')
    print('delta-V: ' + str(round(v_final - v_initial, 1)) + 'm/s')
    print('e_final: ' + str(orbit.ecc.value))

    # df.to_csv('orbit.txt')

    x = df['x'].to_numpy()
    y = df['y'].to_numpy()
    z = df['z'].to_numpy()

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z)

    v = []
    for i in range(len(df)):
        v.append( gnc.magnitude([df['x'].iloc[i], df['y'].iloc[i], df['z'].iloc[i]]) )

    limit = max(v)

    ax.set_xlim([-limit-500, limit+500])
    ax.set_ylim([-limit-500, limit+500])
    ax.set_zlim([-limit-500, limit+500])

    plt.show()
    

    plt.plot(df['time'], df['altitude'])
    plt.grid()
    plt.xlabel('time (sec)')
    plt.ylabel('altitude (km)')
    plt.show()