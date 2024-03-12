# from orbit_propagater import fpropagate as fprop
from gnc_library import fpropagate as fprop
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


R_earth = 6378.1

a_initial = (750 + R_earth) << u.km
e_initial = 10**-4 << u.one
i_initial = 97.7 << u.deg
w_initial = 0 << u.deg
O_initial = 0 << u.deg
nu_initial = 0 << u.deg 

orbit = Orbit.from_classical(Earth, a_initial, e_initial, i_initial, O_initial, w_initial, nu_initial)

v_initial = gnc.magnitude( (orbit.v << u.meter/u.second).value )


f_mag = 0.015
mass = 180
A = 2


df = pd.DataFrame()

print('\nstarting simulation...\n')

xdata, ydata, zdata = [], [], []

t = 0
dt = 60*5
t_max = 25000*60 
while t < t_max:

    f_thrust = spiral_in(f_mag, orbit)
    orbit = fprop(orbit, f_thrust, mass, A, t, dt)

    tempdf = pd.DataFrame({'time':[t],
                           'x':[(orbit.r[0] << u.km).value],
                           'y':[(orbit.r[1] << u.km).value],
                           'z':[(orbit.r[2] << u.km).value],
                           'altitude':[(orbit.a << u.km).value - R_earth],
                           'eccentricity':[orbit.ecc.value],
                           'inclination':[(orbit.inc << u.deg).value]})
    df = pd.concat([tempdf, df])

    print("\r" + str(int(100*t/t_max)) + "% complete", end='')
    t += dt

print('\r100% complete', end='')

v_final = gnc.magnitude( (orbit.v << u.meter/u.second).value )
delta_v = v_final - v_initial

print('\r' + str(int(abs(delta_v))) + ' m/s of delta-V consumed for manuever:')
print('\tspiral in from ' + str(int(df['altitude'].iloc[-1]+1)) + 'km to ' +
      str(int(df['altitude'].iloc[0]+1)) + 'km')
print('\teccentricity wound up being ' + str(round(df['eccentricity'].iloc[0], 4)))
print('\tinclination ended up at ' + str(round(df['inclination'].iloc[0], 2)) + ' degrees')
print('\n')

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

plt.figure()
plt.title('altitude vs. time')
plt.plot(df['time'], df['altitude'])
plt.grid()
plt.xlabel('time (sec)')
plt.ylabel('altitude (km)')

plt.figure()
plt.title('eccentricity vs. time')
plt.plot(df['time'], df['eccentricity'])
plt.grid()
plt.xlabel('time (sec)')
plt.ylabel('eccentricity')

plt.figure()
plt.title('inclination vs. time')
plt.plot(df['time'], df['inclination'])
plt.grid()
plt.xlabel('time (sec)')
plt.ylabel('inclination (deg)')

plt.show()