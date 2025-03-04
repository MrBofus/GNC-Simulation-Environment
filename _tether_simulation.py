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

# import flight software
import flight_software.tether_flight_software.flight_software as fs

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
argp = 0.001    # degrees
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
w_initial = 0.2*np.array([0.016, 0.08, -0.016])

# build the state handler for the simulation
state = gnc.satelliteState(satellite_orbit, moment_of_inertia, satellite_mass,
                           satellite_area, q_initial, w_initial)


# timestep to run physics simulation
physicsHz = 20             # Hz

# timestep to run flight software
flightSoftwareHz = 5      # Hz

# total simulation time
simTime = 30*24*3600   # s


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# reaction wheel parameters
#
# https://www.rocketlabusa.com/assets/Uploads/10-mNms-RW-0.01-Data-Sheet.pdf
# power consumption at max speed / max torque: 1.05 W
# power consumption at max speed / no torque: 0.25 W
# power consumption at no speed / no torque: 0.11 W

max_wheel_speed = 6500*0.1047198        # rad/s
wheel_moi = 0.018 / max_wheel_speed     # Nms^2 or kgm^2
max_wheel_torque = 1*10**-3             # Nm
min_wheel_speed = 0.1*0.1047198         # rad/s

power_consumption_idle = 0.11        # W
power_consumption_active = 1.05      # W

reactionWheelAssembly = gnc.reactionWheelAssembly(wheel_moi, max_wheel_speed, max_wheel_torque, min_wheel_speed)


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# magnetorquer parameters
#
# Jagsat-specific hardware, maybe link to ICD?

n_turns = 2*40                  # unitless (times two to reflect 2 magnetorquers)
length = 20 * 10**-2            # meters
width = 10 * 10**-2             # meters
max_current = 40 * 10 ** -3     # Amps

power_consumption_idle = 0               # W
power_consumption_active = 5*max_current # W

magnetorquerAssembly = gnc.magnetorquerAssembly(n_turns, length, width, max_current)


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# thruster parameters
#
# https://www.enpulsion.com/order/enpulsion-nano/

thrust_magnitude = 330 * 10**-6 # N
propellant_mass = 220 * 10**-3  # kg

power_consumption_idle = 8      # W
power_consumption_active = 40   # W


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# initialize simulation variables

controlTorque = [0, 0, 0]       # initial control torque
disturbanceTorque = [0, 0, 0]   # initial disturbance torque
thrust = [0, 0, 0]              # initial thrust

# initialize dataframe to store results
df = pd.DataFrame()

# initialize counter and time to zero
counter = 0
t = 0 


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# flight software begin {#467, 5}
# begin flight software and update user variables if any

scheduler = fs.schedulerApp()
scheduler._update_user_variables(bDot_gain=10**5, smc_kp=0.003, smc_kd=0.015, gs_range=0,
                                 rg_noise=0, rg_bias=0, st_noise=0, mg_noise=0, gps_noise_r=0, gps_noise_v=0)


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# begin simulation

print('\n')
while t < simTime:
    # add random disturbance torques
    disturbanceTorque = (10**-7)*np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)]) # Nm
    # add random disturbance forces
    disturbanceForce = (10**-5)*np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)]) # N


    # propagate attitude one timestep
    state.propagateAttitude(controlTorque, disturbanceTorque, (1/physicsHz))
    # propagate posiition one timestep
    state.propagatePosition(thrust, disturbanceForce, (1/physicsHz))
    

    # if this timestep is a software timestep, iterate flight software once
    if counter % int(physicsHz / flightSoftwareHz) == 0:

        # advance flight software one timestep
        scheduler.iterate(state)

        # write commands to reaction wheels
        reactionWheelAssembly.commandReactionWheels(scheduler.rw_command)
        # write commands to magnetorquers
        magnetorquerAssembly.commandMangetorquers(scheduler.mt_command)
    

    # allow reaction wheels to influence physical state
    reactionWheelAssembly.actuateReactionWheels(1/physicsHz)
    # allow magnetorquers to influence physical state
    magnetorquerAssembly.actuateMagnetorquers(state)
    

    # consolidate physical influences
    controlTorque = copy.deepcopy(reactionWheelAssembly.wheel_torques) + copy.deepcopy(magnetorquerAssembly.torque)
    # determine thrust
    if scheduler.thruster_command > 0:
        thrust = scheduler.thruster_command * np.array([scheduler.setpoint[0], scheduler.setpoint[1], scheduler.setpoint[2]]) # qm.quaternion_to_axis(state.quaternion)


    # append variables to dataframe to monitor results
    df = gnc.appendDataFrame(df, state, t, scheduler, scheduler.q_error, reactionWheelAssembly)
    
    # update user about simulation status
    if counter%1000 == 0: print('\r', str(int(100*t/simTime)) + '% complete -- (w is ' + str(round(qm.magnitude(state.angularRate), 4)) + 
                                ', a is ' + str(state.orbit.a.value) + ')        ' , end='')

    # advance counter and simulation time
    counter += 1
    t += (1/physicsHz)

    # if target state is reached early, go ahead and end simulation
    if scheduler.mode == 'exit':
        break


# when simulation is finished, stop the timer
t_to_finish = time.monotonic()

# display time to finish simulation
print('\nsimulation took ' + str(round((t_to_finish-t_to_start) / 60, 1)) + ' minutes')

# generate plots for user
gnc.plot_quaternion_error(df)
# show results of simulation
gnc.plot_orbit_transfer(df)

# write results to csv for future analysis
df.to_csv('_out/simulation_results_tether.txt')