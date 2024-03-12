import gnc_core.gnc_library as gnc
import numpy as np
import pandas as pd
import quaternionMath as qm
import time
import random
import copy

import flight_software as fs


t_to_start = time.monotonic()



a = 565
e = 0.001
i = 97.5
raan = 90
argp = 0
nu = -179.9

Ixx = 10
Iyy = 10
Izz = 10

satellite_area = 1.1



satellite_orbit = [(a+6378.1)*10**3, e, i, raan, argp, nu]
moment_of_inertia = [Ixx, Iyy, Izz]
satellite_mass = 25
q_initial = qm.normalize(np.array([random.random(), random.random(), random.random(), random.random()]))
w_initial = np.array([0.06, 0.01, -0.06])


state = gnc.satelliteState(satellite_orbit, moment_of_inertia, satellite_mass,
                           satellite_area, q_initial, w_initial)




# https://nanoavionics.com/cubesat-components/cubesat-reaction-wheels-control-system-satbus-4rw/
max_wheel_speed = 6500*0.1047198
wheel_moi = (20*10**-3) / max_wheel_speed
max_wheel_torque = 3.2*10**-3
min_wheel_speed = 0.1*0.1047198
reactionWheelAssembly = gnc.reactionWheelAssembly(wheel_moi, max_wheel_speed, max_wheel_torque, min_wheel_speed)

n_turns = 40
length = 20 * 10**-2
width = 10 * 10**-2
max_current = 20 * 10 ** -3
magnetorquerAssembly = gnc.magnetorquerAssembly(n_turns, length, width, max_current)


measurements = fs.measured_state()


thrust_magnitude = 5 * 10**1

ctl = [0, 0, 0]
controlTorque = [0, 0, 0]
disturbanceTorque = [0, 0, 0]
thrust = [0, 0, 0]
disturbanceForce = 0


df = pd.DataFrame()


physicsHz = 20            # Hz
flightSoftwareHz = 5      # Hz

counter = 0
t = 0                                     # s
simTime = 0.3*state.orbit.period.value   # s

ff = 0
current_ = 'detumble'

print('\n')
while t < simTime:
    
    disturbanceTorque = (10**-9)*np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)])
    state.propagateAttitude(controlTorque, disturbanceTorque, (1/physicsHz))
    state.propagatePosition(thrust, disturbanceForce, (1/physicsHz))
    
    '''
    v_direction = np.array(gnc.normalize( (state.orbit.v << u.meter / u.second).value ))
    r_direction = np.array(gnc.normalize( (state.orbit.r << u.meter).value ))
    
    limit = 500
    if ((state.orbit.r[2] << u.km).value > - limit and (state.orbit.r[2] << u.km).value < + limit) and ff == 0:
        h_direction = np.cross(v_direction, r_direction)
        unitvec = gnc.normalize( 0.0*v_direction + 1.0*h_direction )
        thrust = thrust_magnitude * np.array( unitvec )
        
    elif ((state.orbit.r[2] << u.km).value > - limit and (state.orbit.r[2] << u.km).value < + limit) and ff == 1:
        h_direction = -np.cross(v_direction, r_direction)
        unitvec = gnc.normalize( 0.0*v_direction + 1.0*h_direction )
        thrust = thrust_magnitude * np.array( unitvec )
    
    elif (state.orbit.r[2] << u.km).value > + limit or (state.orbit.r[2] << u.km).value < - limit:
        unitvec = gnc.normalize( 1.0*v_direction + 0.0*h_direction )
        thrust = (thrust_magnitude/50) * np.array( unitvec )
    
    else:
        thrust = [0, 0, 0]
        
    
    if (state.orbit.r[2] << u.km).value > + limit:
        ff = 0
    elif (state.orbit.r[2] << u.km).value < - limit:
        ff = 1
    '''
    
    
    #############
    

    if counter % int(physicsHz / flightSoftwareHz) == 0:


        #############


        fs._navigation(state, measurements)


        #############

        if current_ == 'detumble':
            mt_command, rw_command, \
                current_, q_error = fs._control_detumble(measurements)
            
        elif current_ == 'dataCollection':
            mt_command, rw_command, \
                current_, q_error = fs._control_dataCollection(measurements)


        reactionWheelAssembly.commandReactionWheels(rw_command)
        magnetorquerAssembly.commandMangetorquers(mt_command)
    

    #############
    
    reactionWheelAssembly.actuateReactionWheels(1/physicsHz)
    magnetorquerAssembly.actuateMagnetorquers(state)
    
    controlTorque = copy.deepcopy(reactionWheelAssembly.wheel_torques) + copy.deepcopy(magnetorquerAssembly.torque)
    # print(controlTorque)
    df = gnc.appendDataFrame(df, state, t, q_error, reactionWheelAssembly)
    
    if counter%1000 == 0:
        print('\r', str(int(100*t/simTime)) + '% complete -- (w is ' + str(round(qm.magnitude(state.angularRate), 4)) + ')        ' , end='')
    
    counter += 1
    t += (1/physicsHz)



t_to_finish = time.monotonic()

# gnc.plot_orbit(df)
gnc.plot_quaternion_error(df)

print('\nsimulation took ' + str(round((t_to_finish-t_to_start) / 60, 1)) + ' minutes')

df.to_csv('__out/simulation_results.txt')