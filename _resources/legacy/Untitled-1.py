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
vvec = -1* np.array([ df['vx'].iloc[0], df['vy'].iloc[0], df['vz'].iloc[0] ])

unitvec = normalize( rvec )

az = np.arctan2(unitvec[0], -unitvec[2])
el = np.arcsin( unitvec[1] )

print(az*180/np.pi)
print(el*180/np.pi)

zvec = np.array([0, 0, 1])
angle = np.arccos(np.dot(zvec, vvec)/(np.sqrt(zvec[0]**2 + zvec[1]**2 + zvec[2]**2)*np.sqrt(vvec[0]**2 + vvec[1]**2 + vvec[2]**2)))
print(angle*180/np.pi)


phi = 10 * 4.8481*10**-6
theta = 10 * 4.8481*10**-6
psi = 100 * 4.8481*10**-6

q1 = np.cos(phi/2)*np.cos(theta/2)*np.cos(psi/2) + np.sin(phi/2)*np.sin(theta/2)*np.sin(psi/2)
q2 = np.sin(phi/2)*np.cos(theta/2)*np.cos(psi/2) - np.cos(phi/2)*np.sin(theta/2)*np.sin(psi/2)
q3 = np.cos(phi/2)*np.sin(theta/2)*np.cos(psi/2) + np.sin(phi/2)*np.cos(theta/2)*np.sin(psi/2)
q4 = np.cos(phi/2)*np.cos(theta/2)*np.sin(psi/2) - np.sin(phi/2)*np.sin(theta/2)*np.cos(psi/2)

print(q1, q2, q3, q4)

import matplotlib.pyplot as plt

plt.plot(df['time'], df['wx'], label='wx')
plt.plot(df['time'], df['wy'], label='wy')
plt.plot(df['time'], df['wz'], label='wz')

plt.grid()
plt.show()

'''
import sys
sys.path.append("../GNC-Simulation-Environment")
import gnc_core.gnc_library as gnc

r = 1000*np.array([ df['x'].iloc[0], df['y'].iloc[0], df['z'].iloc[0] ])

print(r)
lon, lat, alt = gnc.ECI_to_ECEF(r, 100)
print(lon, lat, alt)
r_ = gnc.ECEF_to_ECI(lat, lon, alt, 100)
print(r_)
'''