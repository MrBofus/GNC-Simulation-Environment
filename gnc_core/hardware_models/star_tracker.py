import copy
import random

def pull_star_tracker(state, **kwargs):
    measured_quaternion = copy.deepcopy(state.quaternion)
    for i in range(len(measured_quaternion)):
        measured_quaternion[i] += 2*kwargs['noise']*(0.5 - random.random())
    
    return measured_quaternion