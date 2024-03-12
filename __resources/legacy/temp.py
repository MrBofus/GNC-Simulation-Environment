from orbit_propagater import fpropagate as fprop
# from gnc_library import fpropagate as fprop
import gnc_library as gnc
from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Earth
from astropy import units as u
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def spiral_in(f_mag, orbit):
    altitude = (orbit.a << u.km).value - 6378.1
    if altitude > 565:
        thrust = -1 * f_mag * np.array( gnc.normalize(orbit.v.value) )
    
    else:
        thrust = [0, 0, 0]
    
    return thrust

def spiral_out(f_mag, orbit):
    altitude = (orbit.a << u.km).value - 6378.1
    if altitude < 567:
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

f_mag = 0.001
mass = 80
A = 2


df = pd.DataFrame()

print('\n')

xdata, ydata, zdata = [], [], []

t = 0
dt = 60*1
t_max = 2500*60 
while t < t_max:

    # f_thrust = spiral_out(f_mag, orbit)
    f_thrust = 0
    orbit = fprop(orbit, f_thrust, mass, A, t, dt)

    tempdf = pd.DataFrame({'time':[t],
                           'x':[(orbit.r[0] << u.km).value],
                           'y':[(orbit.r[1] << u.km).value],
                           'z':[(orbit.r[2] << u.km).value],
                           'altitude':[(orbit.a << u.km).value - R_earth],
                           'inclination':[(orbit.inc << u.deg).value]})
    df = pd.concat([tempdf, df])

    print("\r" + str(int(100*t/t_max)) + "% complete", end='')
    t += dt

print('\n')
v_final = gnc.magnitude((orbit.v << u.meter/u.second).value)

# print(v_final - v_initial)
print('v_initial: ' + str(round(v_initial, 1)) + 'm/s')
print('v_final: ' + str(round(v_final, 1)) + 'm/s')
print('delta-V: ' + str(round(v_final - v_initial, 1)) + 'm/s')
print('e_final: ' + str(orbit.ecc.value))

df.to_csv('orbit.txt')

'''
x = np.array(df['x'])
y = np.array(df['y'])
z = np.array(df['z'])

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x, y, z)

v = []
for i in range(len(df)):
    v.append( gnc.magnitude([df['x'].iloc[i], df['y'].iloc[i], df['z'].iloc[i]]) )

limit = max(v)

ax.set_xlim([-limit-5000, limit+5000])
ax.set_ylim([-limit-5000, limit+5000])
ax.set_zlim([-limit-5000, limit+5000])

plt.show()
'''

plt.plot(df['time'], df['altitude'])
plt.grid()
plt.xlabel('time (sec)')
plt.ylabel('altitude (km)')
plt.show()