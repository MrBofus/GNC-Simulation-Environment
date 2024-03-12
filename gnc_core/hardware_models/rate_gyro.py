import copy
import random

def pull_gyro(state, **kwargs):
    measured_angular_rate = copy.deepcopy(state.angularRate)
    for i in range(len(measured_angular_rate)):
        measured_angular_rate[i] += 2*kwargs['noise']*(0.5 - random.random())
    
    return measured_angular_rate