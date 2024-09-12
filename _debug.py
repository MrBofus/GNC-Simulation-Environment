# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# written by: ME !! :)
#
# import libraries

import gnc_core.gnc_library as gnc
import numpy as np
import pandas as pd
import quaternion_math.quaternionMath as qm
import time
import random
import copy


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


# satellite dimensions ( see above )
satellite_x = 20 * 10**-2  # m
satellite_y = 20 * 10**-2  # m
satellite_z = 20 * 10**-2 # m

# satellite mass
satellite_mass = 1.0 # kg

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
w_initial = 0.001*np.array([0.016, 0.08, -0.016])

# build the state handler for the simulation
state = gnc.satelliteState(satellite_orbit, moment_of_inertia, satellite_mass,
                           satellite_area, q_initial, w_initial)


# timestep to run physics simulation
physicsHz = 5             # Hz

# timestep to run flight software
flightSoftwareHz = 1      # Hz

# total simulation time
simTime = 0.2*state.orbit.period.value   # s


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# magnetorquer parameters

n_turns = 20                    # unitless
length = 20 * 10**-2            # meters
width = 20 * 10**-2             # meters
max_current = 40 * 10 ** -3     # Amps

power_consumption_idle = 0               # W
power_consumption_active = 5*max_current # W

magnetorquerAssembly = gnc.magnetorquerAssembly(n_turns, length, width, max_current)


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


integrated_error = np.array([0., 0., 0., 0.])

# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# begin simulation

print('\n')
while t < simTime:
    # add random disturbance torques
    disturbanceTorque = 0*np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)]) # N

    # propagate attitude one timestep
    state.propagateAttitude(controlTorque, disturbanceTorque, (1/physicsHz))

    # propagate posiition one timestep
    state.propagatePosition(thrust, disturbanceForce, (1/physicsHz))

    

    # if this timestep is a software timestep, iterate flight software once
    if counter % int(physicsHz / flightSoftwareHz) == 0:

        '''
        # kp = [0.03, 0.03, 0.03]
        # kp = 0.03
        # kd = 5

        # kp = 0.003
        # ki = 0.00018
        # kd = 50
        
        tx = 10
        
        kx = ((4*np.pi**2)*Ixx)/(tx**2) - (Izz - Iyy)*(0)**2
        ky = ((4*np.pi**2)*Iyy)/(tx**2) - 4*(Izz - Ixx)*(0)**2
        kz = ((4*np.pi**2)*Izz)/(tx**2) - 3*(Iyy - Ixx)*(0)**2
        
        kdx = 16*np.pi*1 * Ixx / tx
        kdy = 16*np.pi*1 * Iyy / tx
        kdz = 16*np.pi*1 * Izz / tx
        
        kp = [kx, ky, kz]
        kd = [kdx, kdy, kdz]
        
        print(kp, kd)
        
        ki = 0
        '''

        kp = 0.03
        kd = 5
        
        quaternion = state.quaternion
        setpoint = np.array([0, 0, 0, 1])
        q_err = qm.quaternionDifference(quaternion, setpoint)
        
        rate_error = np.array(state.angularRate) - np.array([0, 0, 0])

        required_torque_for_control = []
        for i in range(3):
            
            if abs(q_err[i]) < 0.1:
                u = 0.1*kp * q_err[i] + 10*kd * rate_error[i]
            else:
                u = kp * q_err[i] + kd * rate_error[i]
            
            required_torque_for_control.append(u)
        
        current_command = qm.return_c(state.B_body, np.array(required_torque_for_control))

        # write commands to magnetorquers
        magnetorquerAssembly.commandMangetorquers(current_command)
    

    
    # allow magnetorquers to influence physical state
    magnetorquerAssembly.actuateMagnetorquers(state)
    
    # consolidate physical influences
    controlTorque = magnetorquerAssembly.torque


    # append variables to dataframe to monitor results
    df = gnc.appendDataFrame_spec(df, state, t, q_err)
    
    # update user about simulation status
    if counter%1000 == 0:
        print('\r', str(int(100*t/simTime)) + '% complete -- (w is ' + str(round(qm.magnitude(state.angularRate), 4)) + ')        ' , end='')
    
    # advance counter and simulation time
    counter += 1
    t += (1/physicsHz)


# when simulation is finished, stop the timer
t_to_finish = time.monotonic()

# generate plots for user
gnc.plot_spec(df)

# display time to finish simulation
print('\nsimulation took ' + str(round((t_to_finish-t_to_start) / 60, 1)) + ' minutes')

# write results to csv for future analysis
df.to_csv('_out/simulation_results_spec.txt')