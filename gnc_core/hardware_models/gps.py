import copy
import random
from astropy import units as u

def pull_gps(state, **kwargs):
    rvec = copy.deepcopy((state.orbit.r << u.km).value)
    vvec = copy.deepcopy((state.orbit.v << u.km/u.second).value)
    
    for i in range(len(rvec)):
        rvec[i] += 2*kwargs['r_noise']*(0.5 - random.random())
        vvec[i] += 2*kwargs['v_noise']*(0.5 - random.random())
    
    return rvec, vvec