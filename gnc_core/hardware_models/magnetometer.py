import copy
import random

def pull_magnetometer(state, **kwargs):
    measured_Bfield = copy.deepcopy(state.B_body)
    for i in range(len(measured_Bfield)):
        measured_Bfield[i] += 2*kwargs['noise']*(0.5 - random.random())
    
    return measured_Bfield