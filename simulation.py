import gnc_library as gnc
import numpy as np
import pandas as pd
import time
from astropy import units as u

t_to_start = time.monotonic()



a = 565
e = 0.001
i = 51.6 * 0 + 10
raan = 90
argp = 0
nu = 0

Ixx = 10
Iyy = 10
Izz = 10

satellite_area = 1.1



satellite_orbit = [(a+6378.1)*10**3, e, i, raan, argp, nu]
moment_of_inertia = [Ixx, Iyy, Izz]
satellite_mass = 25
q_initial = gnc.normalize(np.array([1, 2, 3, 1]))
w_initial = np.array([0.001, 0.001, 0.001])


state = gnc.satelliteState(satellite_orbit, moment_of_inertia, satellite_mass,
                           satellite_area,  q_initial, w_initial)



kp = 1
kd = 10
sigma = 0.05
order = 3

p2Controller = gnc.p2Controller(kp, kd)
smcController = gnc.smcController(kp, kd, sigma, order)

max_wheel_speed = 100
wheel_moi = 0.01 / max_wheel_speed
max_wheel_torque = 0.01
min_wheel_torque = 0
reactionWheelAssembly = gnc.reactionWheelAssembly(wheel_moi, max_wheel_speed, max_wheel_torque, min_wheel_torque)



thrust_magnitude = 5 * 10**1

ctl = [0, 0, 0]
controlTorque = [0, 0, 0]
disturbanceTorque = [0, 0, 0]
thrust = [0, 0, 0]
disturbanceForce = 0


df = pd.DataFrame()


physicsHz = 50            # Hz
flightSoftwareHz = 5        # Hz

counter = 0
t = 0                                    # s
simTime = 0.35*state.orbit.period.value    # s

ff = 0

while t < simTime:
    
    disturbanceTorque = (10**-9)*np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)])
    state.propagateAttitude(controlTorque, disturbanceTorque, (1/physicsHz))
    # state.propagatePosition(thrust, disturbanceForce, (1/physicsHz))
    
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
    
    
    
    
    if counter % int(physicsHz / flightSoftwareHz) == 0:
        
        p2Controller.input_function(state.angularRate, state.quaternion)
        # smcController.input_function(state.angularRate, state.quaternion)
        
        ctl = p2Controller.output_function( [0, 0, 0, 1] )
        # ctl = smcController.output_function( [0, 0, 0, 1] )
    
    
    
    reactionWheelAssembly.actuateReactionWheels(ctl, (1/physicsHz))
    
    
    controlTorque = reactionWheelAssembly.wheel_torques
    
    df = gnc.appendDataFrame(df, state, t, p2Controller.quaternion_error)
    
    if counter%1000 == 0:
        print(str(int(100*t/simTime)) + '% complete')
    
    counter += 1
    t += (1/physicsHz)

t_to_finish = time.monotonic()

# gnc.plot_orbit(df)
gnc.plot_quaternion_error(df)

print('simulation took ' + str(round((t_to_finish-t_to_start) / 60, 1)) + ' minutes')