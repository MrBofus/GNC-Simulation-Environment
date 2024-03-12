from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Earth
from astropy import units as u
# from pyatmos import download_sw_jb2008, read_sw_jb2008
# from pyatmos import jb2008
import gnc_library as gnc
import numpy as np


# swfile = download_sw_jb2008() 

'''
def compute_drag(r, v, A, m, t):
    swdata = read_sw_jb2008(swfile)
    
    time = '2014-07-22 22:18:45'
    lon, lat, alt = gnc.ECI_to_ECEF(r, t)
    jb08 = jb2008(time, (lat, lon, alt/10**3), swdata)
    
    rho = jb08.rho
    Cd = 2.2
    
    return (1/2) * rho * Cd * A * v**2
'''

def integrate_accel(accel, velocity, dt):
    v_candidate = []
    for i in range(3):
        v_candidate.append(velocity[i] + accel[i]*dt)
    
    return v_candidate



def fpropagate(orbit, f_thrust, mass, A, t, dt):
    
    r, v = (orbit.r << u.meter).value, (orbit.v << u.meter/u.second).value
    
    # f_drag_magnitude = compute_drag(r, v, A, mass, t)
    # f_drag = -1 * f_drag_magnitude * np.array(gnc.normalize(v))
    
    
    # f_total = f_drag + f_thrust
    f_total = np.array(f_thrust)
    
    v = integrate_accel(f_total/mass, v, dt)
    
    return Orbit.from_vectors(Earth, r << u.meter, v << u.meter/u.second)