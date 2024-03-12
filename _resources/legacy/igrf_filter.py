##````````````````````````````````````````````````````````````````````````````````````````````````````````##
# import libraries

# standard imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# goofy libraries
import pyIGRF as igrf
from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Earth
from astropy import units as u
import pyproj

##````````````````````````````````````````````````````````````````````````````````````````````````````````##
# define some helper functions

# the IGRF model takes lat/long as inputs, so we need a
# way of converting the position of the satellite to
# lat_long
def ECI_to_LatLong(r):
    
    # account for rotation of the Earth
    gamma = (360/(23.9345*3600))
    gamma *= np.pi/180
    
    # compute x, y, and z position of satellite over Earth
    x = r[0]*np.cos(gamma) - r[1]*np.sin(gamma) # meters
    y = r[0]*np.sin(gamma) + r[1]*np.cos(gamma) # meters
    z = r[2] # meters
    
    # build pyproj transformer based on WGS84 Ellipsoid
    transformer = pyproj.Transformer.from_crs({"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
                                              {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'})
    
    # compute lat/lon/altitude in degrees/km
    lon, lat, alt = transformer.transform(x, y, z, radians=False)
    
    return lon, lat, alt


# given satellite position, determine the true magnetic field
# given the IGRF model.
def return_true_magnetic_field(position):

    # determine lat/lon/alt given position
    lon, lat, alt = ECI_to_LatLong(position)

    # compute igrf value given position and year
    v = igrf.igrf_value(lat, lon, alt/1000, 2022)

    # convert from nT to T and return magnetic field
    return (10**-9)*np.array([v[3], v[4], v[5]])


# given satellite position, measure the magnetic field.
# adds sensor noise to show error-filled measurement.
def measure_magnetic_field(position):

    # determine true magnetic field
    true_field = return_true_magnetic_field(position)

    # add random sensor noise
    measured_field = np.array([ true_field[0] + np.random.normal(0, 10**-6),
                                true_field[1] + np.random.normal(0, 10**-6),
                                true_field[2] + np.random.normal(0, 10**-6) ])
    
    # return measurement
    return measured_field

##````````````````````````````````````````````````````````````````````````````````````````````````````````##
# define the filter

def alpha_beta_gamma_filter(measurement, position, weight):
    
    #````````````````````````````````````````````````````````#
    # Predict state using the IGRF equation:

    # determine lat/lon/alt given position
    lon, lat, alt = ECI_to_LatLong(position)

    # compute igrf value given position and year
    v = igrf.igrf_value(lat, lon, alt/1000, 2022)

    # convert from nT to T and package into one vector
    x_n = (10**-9)*np.array([v[3], v[4], v[5]])

    # we will use x_n as our predicted state

    #````````````````````````````````````````````````````````#

    # putting our measurement in the same formalism:
    z_n = measurement

    # we now have predicted state from the theoretical approach
    # and our measured state from our sensor.

    #````````````````````````````````````````````````````````#
    # the filter:

    x_N = x_n + weight*(z_n - x_n)

    # where x_N is our estimated state given the prediction
    # and the measurement.

    # ** as far as I know we don't need the explicit
    #    gamma term as that's covered inside the
    #    igrf model? functionally speaking

    return x_N


##````````````````````````````````````````````````````````````````````````````````````````````````````````##
# initialize variables

# define radius of the Earth for later
R_earth = 6378.1 # km

a_initial = (565 + R_earth) << u.km     # semimajor axis in km
e_initial = 10**-4 << u.one             # eccentriciy
i_initial = 45 << u.deg                 # inclination in degrees
w_initial = 0 << u.deg                  # argument of periapsis in degrees
O_initial = 0 << u.deg                  # RAAN in degrees
nu_initial = 0 << u.deg                 # true anomaly in degrees

# build the orbit object of the satellite
orbit = Orbit.from_classical(Earth, a_initial, e_initial, i_initial, O_initial, w_initial, nu_initial)

dt = 10         # define timestep in seconds
t_max = 5400    # simulation runs for 5400 seconds, or one orbit
t = 0           # initialize time to zero

counter = 0     # initialze counter to zero

df = pd.DataFrame() # initialize dataframe

print('\n')
##````````````````````````````````````````````````````````````````````````````````````````````````````````##
# begin simulation
while t < t_max:

    orbit = orbit.propagate(dt << u.second) # propagate orbit by 1 timestep

    # make a meaurement of position using gps:
    position = (orbit.r << u.meter).value # (in meters)

    # measure magnetic field with magnetometer
    measured_field = measure_magnetic_field(position)

    # filter this measurement with filter:
    # vary the weight-- closer to one means trust the IGRF model more,
    #                   closer to 0 means trust the sensor more.
    filtered_field = alpha_beta_gamma_filter(measured_field, position, 0.3)

    # append measurements to dataframe
    tempdf = pd.DataFrame({'time':[t], 
                           'measured field':[measured_field[0]],
                           'filtered field':[filtered_field[0]]})
    df = pd.concat([df, tempdf])

    # every so often print update:
    if counter % 100 == 0: print("\r" + str(int(100*t/t_max)) + "% complete", end='')

    # increment time
    t += dt
    counter += 1

print("\r100% complete")

# make graph showing measured and filtered data
plt.title('measured field and filtered field vs. time')
plt.plot(df['time'], df['measured field'], label='measured')
plt.plot(df['time'], df['filtered field'], label='filtered')
plt.grid()
plt.show()