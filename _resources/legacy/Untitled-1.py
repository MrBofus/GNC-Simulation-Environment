# import sys
# print(sys._getframe().f_back.f_code.co_name)

'''
class measured_state():
    def __init__(self):
        self.m = {}

    def __repr__(self, key):
        return self.m[0]

    def update(self, **kwargs):
        self.m = {}
        for arg, value in kwargs.items():
            self.m.update({arg:value})


m = measured_state()

m.update(a=1, b=2, c=3)

print(m.m['a'])
print(m.m['b'])
print(m.m['c'])

m.update(a=3, b=4, c=5)

print(m.m['a'])
print(m.m['b'])
print(m.m['c'])
'''

import pandas as pd
import numpy as np

def normalize(vec):
    m = 0
    for v in vec:
        m += v**2
    
    v_ = []
    for v in vec:
        v_.append(v/np.sqrt(m))
    
    return np.array(v_)

df = pd.read_csv('_out/simulation_results.txt')

q = np.array([ df['q1'].iloc[0], df['q2'].iloc[0], df['q3'].iloc[0], df['q4'].iloc[0] ])
rvec = -1* np.array([ df['x'].iloc[0], df['y'].iloc[0], df['z'].iloc[0] ])

unitvec = normalize( rvec )

az = np.arctan2(unitvec[0], -unitvec[2])
el = np.arcsin( unitvec[1] )

print(az*180/np.pi)
print(el*180/np.pi)