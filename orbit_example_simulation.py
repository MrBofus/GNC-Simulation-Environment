# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# written by: ME !! :)
#
# import libraries

import gnc_core.gnc_library as gnc
import numpy as np
import pandas as pd
import time

# import flight software
import example_flight_software.example_flight_software_orbit.flight_software as fs

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



# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# simulation parameters {#f91, 22}
# simulation parameters

# build the state handler for the simulation (use special physics object for orbit-only simulation)
state = gnc.satelliteState_orbit(satellite_orbit, satellite_mass, satellite_area)


# timestep to run physics simulation
physics_timestep = 3*60                 # seconds


# maximum simulation time (this example will cut off early once target state is reached)
simTime = 100*24*3600   # s


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

scheduler = fs.schedulerApp()
scheduler._update_user_variables(f_mag = 330 * 10**-6)


# `````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````` #
# begin simulation

print('\n')
while t < simTime:
    # add random disturbance forces
    disturbanceForce = (10**-5)*np.array([np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)]) # N

    # propagate posiition one timestep
    state.propagatePosition(thrust, disturbanceForce, physics_timestep)
    


    # iterate flight software one timestep
    scheduler._iterate(state)
    
    # pass command to thruster from flight software
    thrust = scheduler.thruster_command



    # append variables to dataframe to monitor results
    df = gnc.appendDataFrame_orbit(df, state, scheduler, t)
    
    # update user about simulation status
    if counter%1000 == 0:
        print('\r', str(int(100*t/simTime)) + '% complete         ' , end='')
    
    # advance counter and simulation time
    counter += 1
    t += physics_timestep

    # if target state is reached early, go ahead and end simulation
    if scheduler.mode == 'exit':
        break


# when simulation is finished, stop the timer
t_to_finish = time.monotonic()

# display time to finish simulation
print('\nsimulation took ' + str(round((t_to_finish-t_to_start) / 60, 1)) + ' minutes')

# show results of simulation
gnc.plot_orbit_transfer(df)

# write results to csv for future analysis
df.to_csv('_out/simulation_results.txt')