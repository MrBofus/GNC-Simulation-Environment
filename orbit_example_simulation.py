# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# written by: ME !! :)
#
# import libraries

import gnc_core.gnc_library as gnc
import numpy as np
import pandas as pd
import quaternionMath as qm
import time
import random
from astropy import units as u
import matplotlib.pyplot as plt

# import flight software
import flightsoftware_core.flight_software as fs
from flightsoftware_core.QLaw import solveQLaw

# start simulation timer
t_to_start = time.monotonic()


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# orbit parameters
#
# remember:
#   sun-synchronous: a = 565km, i = 97.5, raan = 157.5, e = 0.0005
#   ISS:             a = 417km, i = 51.6, raan = 77.6,  e = 0.0006

a = 417         # km
e = 0.0006      # unitless
i = 51.6        # degrees
raan = 77.6     # degrees
argp = 0        # degrees
nu = 0          # degrees

# define satellite orbit given parameters
satellite_orbit = [(a+6378.1)*10**3, e, i, raan, argp, nu]


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# satellite physical parameters
#
#                x                               z
#           ___________               _________________________                  
#          \           \             \                         \
#          \           \    y        \                         \    y 
#          \           \             \                         \
#          \___________\             \_________________________\
#
#
#                x
#           ___________          
#          \           \                   
#          \           \                 
#          \           \                 
#          \           \    z             
#          \           \                
#          \           \                                            
#          \           \    
#          \___________\ 
#
#           

# satellite dimensions ( see above )
satellite_x = 10 * 10**-2 # m
satellite_y = 10 * 10**-2 # m
satellite_z = 20 * 10**-2 # m

# satellite mass
satellite_mass = 1.5 # kg

# average area, computed using the min. and max areas
satellite_area = 0.5 * (satellite_x**2 + satellite_z**2)    # m^2 

# moments of inertia ( calculated using uniform assumption )
Ixx = (1/12) * satellite_mass * (satellite_y**2 + satellite_z**2) 
Iyy = (1/12) * satellite_mass * (satellite_x**2 + satellite_z**2)
Izz = (1/12) * satellite_mass * (satellite_x**2 + satellite_y**2)
moment_of_inertia = [Ixx, Iyy, Izz] # kg-m^2


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# simulation parameters {#f91, 22}
# simulation parameters

# initial quaternion ( random )
q_initial = qm.normalize(np.array([random.random(), random.random(), random.random(), random.random()]))

# initial angular rate ( rad/s )
w_initial = 0.2*np.array([0, 0, 0])

# build the state handler for the simulation
state = gnc.satelliteState(satellite_orbit, moment_of_inertia, satellite_mass,
                           satellite_area, q_initial, w_initial)


# timestep to run physics simulation
physicsHz = (1/(3*60))            # Hz


# total simulation time
simTime = 200*state.orbit.period.value   # s


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# thruster parameters
#
# https://www.enpulsion.com/order/enpulsion-nano/

thrust_magnitude = 330 * 10**-6 # N
propellant_mass = 220 * 10**-3 # kg

power_consumption_idle = 8      # W
power_consumption_active = 40   # W


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# initialize simulation variables

controlTorque = [0, 0, 0]       # initial control torque
disturbanceTorque = [0, 0, 0]   # initial disturbance torque
thrust = [0, 0, 0]              # initial thrust
disturbanceForce = 0            # initial disturbance force (unused)

# initialize dataframe to store results
df = pd.DataFrame()

# initialize counter and time to zero
counter = 0
t = 0 


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# flight software begin {#467, 5}
# begin flight software and update user variables if any


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# begin simulation

print('\n')
while t < simTime:
    # add random disturbance torques
    disturbanceTorque = (10**-7)*np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)]) # N

    # propagate posiition one timestep
    state.propagatePosition(thrust, disturbanceForce, (1/physicsHz))

    
    Wp = 1
    rp_min = 400000
    f_mag = thrust_magnitude
    Wa = 10
    We = 10
    Wi = 1
    Wraan = 0
    Wargp = 0
    aT = (530+6378.1)*10**3
    eT = 0.001
    iT = 51.9 * np.pi/180
    raanT = 0
    argpT = 0
    
    gvec, q_dot = solveQLaw(state.orbit, Wp, rp_min, f_mag,
                            Wa, We, Wi, Wraan, Wargp,
                            aT, eT, iT, raanT, argpT)
    
    r_mag = qm.normalize(state.orbit.r.value)
    v_mag = qm.normalize(state.orbit.v.value)
    h_mag = np.cross(v_mag, r_mag)
    
    gvec_local = [r_mag[0]*gvec[0] + v_mag[0]*gvec[1] + h_mag[0]*gvec[2],
                  r_mag[1]*gvec[0] + v_mag[1]*gvec[1] + h_mag[1]*gvec[2],
                  r_mag[2]*gvec[0] + v_mag[2]*gvec[1] + h_mag[2]*gvec[2]]

    '''
    if (state.orbit.nu << u.deg).value > 90 and (state.orbit.nu << u.deg).value < 180:
        thrust = thrust_magnitude * np.array(gvec_local)
    else:
        thrust = [0, 0, 0]
    '''
    thrust = thrust_magnitude * np.array(gvec_local)


    # append variables to dataframe to monitor results
    tempdf = pd.DataFrame({'t':[t], 'a':[(state.orbit.a << u.km).value],
                           'e':[state.orbit.ecc.value], 'i':[(state.orbit.inc << u.deg).value],
                           'q dot':[q_dot]})
    df = pd.concat([df, tempdf])
    
    # update user about simulation status
    if counter%1000 == 0:
        print('\r', str(int(100*t/simTime)) + '% complete -- (w is ' + str(round(qm.magnitude(state.angularRate), 4)) + ')        ' , end='')
    
    # advance counter and simulation time
    counter += 1
    t += (1/physicsHz)


# when simulation is finished, stop the timer
t_to_finish = time.monotonic()

# display time to finish simulation
print('\nsimulation took ' + str(round((t_to_finish-t_to_start) / 60, 1)) + ' minutes')

# generate plots for user
plt.figure(1)
plt.plot([0, df['t'].iloc[-1]], [aT/(10**3)-6371.8, aT/(10**3)-6371.8], 'k-')
plt.plot(df['t'], df['a']-6371.8)
plt.grid()

plt.figure(2)
plt.plot([0, df['t'].iloc[-1]], [eT, eT], 'k-')
plt.plot(df['t'], df['e'])
plt.grid()

plt.figure(3)
plt.plot([0, df['t'].iloc[-1]], [iT*180/np.pi, iT*180/np.pi], 'k-')
plt.plot(df['t'], df['i'])
plt.grid()

plt.figure(4)
plt.plot(df['t'], df['q dot'])
plt.grid()

plt.show()

# write results to csv for future analysis
df.to_csv('_out/simulation_results.txt')