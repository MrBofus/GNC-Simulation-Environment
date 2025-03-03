import gnc_core.hardware_models.magnetometer as mm
import gnc_core.hardware_models.rate_gyro as rg
import gnc_core.hardware_models.star_tracker as st
import gnc_core.hardware_models.gps as gps

from flight_software.tether_flight_software.local_gnc_library import *
from flight_software.tether_flight_software.QLaw import solveQLaw

import copy
import quaternion_math.quaternionMath as qm
from astropy import units as u
import numpy as np


# define printout function specifically for flight software apps,
# that prints in green so the user knows what the flight software
# prints vs. the simulation
def printout(str):
    print('\n\n\033[0;32m  from flight_software.schedulerApp._update_user_variables():\n' + str + '\033[0m')


# this is the actual flight software that runs on the simulated satellite.
class schedulerApp():    
    # define the measured state class, which handles the measurements
    # and puts them in an organized way to access later
    class measured_state():
        def __init__(self):
            self.measurements = {}

        def update(self, **kwargs):
            self.measurements = {}
            for arg, value in kwargs.items():
                self.measurements.update({arg:value})


    # ``````````````````````````````````````````````````````````````````````````````````````````````
    # flight software initialization
    #
    # initialize the flight software with some defaults that we can
    # update later
    def __init__(self, **kwargs):
        if 'timestep' in kwargs: dt=kwargs['timestep']
        else: dt = 3*60
        
        self.defaults = {
                        # ````````````````````````````````````````
                        # defualts for sensors
                        'bDot_gain':0.1, # gain for the b-dot controller for detumble
                
                        'smc_kp':1.0,    # p-gain for the sliding-mode controller, for pointing with wheels
                        'smc_kd':10.0,   # d-gain for the sliding-mode controller, for pointing with wheels
                        'smc_sigma':0.1, # window for sliging mode controller
                        'smc_order':3,   # sliding mode controller order

                        'rg_noise':1/5000, # gaussian sensor noise for rate gyro (rad/s)
                        'rg_bias':0.001,   # bias in rate gyro measurements (rad/s)

                        'st_noise':1/5000, # star tracker sensor noise (quaternion error)
                        
                        'mg_noise':10**-11, # magnetometer sensor noise (T)
                        
                        'gps_noise_r':0.01, # innacuracy in gps position (km)
                        'gps_noise_v':0.01, # innacuracy in gps velocity (km/s)

                        'gs_range':800*10**3, # minimum ground station distance for downlink (m)

                        # ````````````````````````````````````````
                        # defualts for transfer law
                        'Wp':1,             # p-gain for Q-Law
                        'rp_min':(400000 + 6378.1) * 10**3, # minimum periapsis for Q-Law

                        'f_mag':330 * 10**-6, # magnitude of thruster (N)

                        'Wa':10,    # semimajor axis gain for Q-Law
                        'We':10,    # eccentriciy gain for Q-Law
                        'Wi':1,     # inclanation gain for Q-Law
                        'Wraan':0,  # RAAN gain for Q-Law
                        'Wargp':0,  # arguement of periapsis gain for Q-Law

                        'aT':(530 + 6378.1)*10**3, # target semimajor axis (m)
                        'eT':0.001,                # target eccentricity (unitless)
                        'iT':51.9 * np.pi/180,     # target inclination (radian)
                        'raanT':0,                 # target RAAN (radian)
                        'argpT':0,                 # target argument of periapsis (radian)
                        
                        'aT_accuracy':0.0001,      # semimajor axis target threshold (%)
                        'iT_accuracy':0.00001,     # inclination target threshold (%)
                        'eT_accuracy':0.0001,      # eccentriciy target threshold (%)
                        
                        'wakeup_time':2 * (24*3600) / dt,   # time spent in bootup
                        'hold_1_window':2 * (24*3600) / dt, # time spent at hold 1 target
                        'hold_2_window':2 * (24*3600) / dt, # time spent at hold 2 target

                        }
        
        self._p = self.defaults
        self.m = self.measured_state()
        self.systemClock = 0


    # ``````````````````````````````````````````````````````````````````````````````````````````````
    # update user variables
    #
    # function that updates user variables according to the default list.
    # if the updated variable is in the default list, it will update,
    # otherwise it will return a warning and ignore the variable
    def _update_user_variables(self, **kwargs):
        outstr = ''
        for arg, value in kwargs.items():
            if arg in self.defaults:
                outstr += '\t--- updated: ' + arg + '  ' + str(self.defaults[arg]) + ' --> ' + str(value) + '\n'
                self._p[arg]=value
            else:
                outstr += '\033[0;31m\t--- warning: ' + arg + ' not recognized as defualt variable, no action taken\033[0;32m\n'
        
        outstr += '\n\t' + str(self._p)
        printout(outstr)


    # ``````````````````````````````````````````````````````````````````````````````````````````````
    # interface commands
    #
    # write command to simulated reaction wheels
    def _write_command_to_reactionWheel(self, command):
        self.rw_command = command
    
    # write command to simulated magnetorquers
    def _write_command_to_magnetorquer(self, command):
        self.mt_command = command

    # write command to simulated thruster
    def _write_command_to_thruster(self, command):
        self.thruster_command = command

    # write variables to system state
    def _write_variables_to_state(self, mode, q_error):
        if not mode == None:
            self.mode = mode
        if not q_error == None:
            self.q_error = q_error
    

    # ``````````````````````````````````````````````````````````````````````````````````````````````
    # main function
    #
    # main function that executes every physical timestep
    def iterate(self, state):
        # update stored physical state to expose to simulated hardware
        self.state = copy.deepcopy(state)
        # incriment the system clock
        self.systemClock += 1
        
        # initialize in 'init' mode
        if self.systemClock == 1:
            self.mode = 'init'
            self.q_error = [1, 1, 1, 1]
        
        # run guidance, navigation, and control
        self._run_navigation() # determine physical state
        self._run_guidance()   # decide what to do given physical state
        self._run_control()    # do something given the decision from guidance

    # ``````````````````````````````````````````````````````````````````````````````````````````````
    # helper functions
    #
    # check if current orbit is within threshold
    def _validate_orbit(self):
        if (self.state.orbit.a << u.meter).value > (1-self._p['aT_accuracy'])*self._p['aT'] and (self.state.orbit.a << u.meter).value < (1+self._p['aT_accuracy'])*self._p['aT']:
            if (self.state.orbit.inc << u.radian).value > (1-self._p['iT_accuracy'])*self._p['iT'] and (self.state.orbit.inc << u.radian).value < (1+self._p['iT_accuracy'])*self._p['iT']:
                if self.state.orbit.ecc.value > (1-self._p['eT_accuracy'])*self._p['eT'] and self.state.orbit.ecc.value < (1+self._p['eT_accuracy'])*self._p['eT']:
                    return True
                else: return False
            else: return False
        else: return False
        
    # ``````````````````````````````````````````````````````````````````````````````````````````````
    # flight software apps
    #
    # guidance app
    # run guidance to determine what the gnc system should be doing given the current
    # physical state of the spacecraft, as determined by the sensors
    def _run_guidance(self):

        #````````````````````````````
        # boot into 'init' mode, where system waits for a timeout before 
        # proceeding to detumble mode
        if self.mode == 'init':
            if self.systemClock > self._p['wakeup_time']:
                self._write_variables_to_state('detumble', None)
            else:
                pass
        
        #````````````````````````````
        # when the system init timeout is tripped, begin detumble
        elif self.mode == 'detumble':
            # if the magnitude of angular rate is less than 0.002 rad/s,
            # exit detumble mode and enter a parking orbit to decide
            # what to do next
            if qm.magnitude(self.m.measurements['angularRate']) < 0.002:
                printout('detumbled, entering parking orbit...')

                # update the mode
                self._write_variables_to_state('parking_orbit', None)
                # mark the time of detumble
                self.t_crit = self.systemClock
            else:
                pass
        
        #````````````````````````````
        # after detumble, enter a parking orbit
        elif self.mode == 'parking_orbit':
            # after a certain amount of time, enter orbit transfer
            if self.systemClock >  self.t_crit + self._p['wakeup_time']:
                printout('hold complete, entering transfer 1...')

                # update the mode
                self._write_variables_to_state('transfer_1', None)
            else:
                pass

        #````````````````````````````
        # begin first orbit transfer
        elif self.mode == 'transfer_1':
            # check if each orbit parameter is within user-defined threshold
            if self._validate_orbit():
                printout('transfer 1 complete, holding...')

                # update the mode
                self._write_variables_to_state('hold_1', None)
                # mark the time of orbit transfer
                self.t_crit = self.systemClock
            else:
                pass
        
        #````````````````````````````
        # once first target orbit is reached, hold at interim orbit
        elif self.mode == 'hold_1':
            # check how long spacecraft has held in orbit and 
            # once timeout is reached, initiate next transfer
            if self.systemClock > self.t_crit + self._p['hold_1_window']:
                self._write_variables_to_state('exit', None) ##################################
                printout('hold complete, beginning transfer 2...')

                # update user variables to define the next transfer
                self._update_user_variables(aT=self._p['aT'] + 120*10**3, 
                                            iT=self._p['iT'] + 0.1*np.pi/180)
                # update the mode
                self._write_variables_to_state('transfer_2', None)
            else:
                pass
        
        #````````````````````````````
        # initiate the next transfer
        elif self.mode == 'transfer_2':
            # check if each orbit parameter is within user-defined threshold
            if self._validate_orbit():
                printout('transfer 2 complete, holding...')

                # update the mode
                self._write_variables_to_state('hold_2', None)
                # mark the time of orbit transfer
                self.t_crit = self.systemClock
            else:
                pass
        
        #````````````````````````````
        # once second target orbit is reached, hold at interim orbit
        elif self.mode == 'hold_2':
            # check how long spacecraft has held in orbit and 
            # once timeout is reached, initiate next transfer
            if self.systemClock > self.t_crit + self._p['hold_2_window']:
                # exit simulation
                self._write_variables_to_state('exit', None)
            else:
                pass

    # navigation app
    # run navigation to determine the current physical state of the spacecraft
    # and return useful information depending on the mode
    def _run_navigation(self):
        wvec = rg.pull_gyro(self.state, noise=self._p['rg_noise'])
        qvec = st.pull_star_tracker(self.state, noise=self._p['st_noise'])
        bvec = mm.pull_magnetometer(self.state, noise=self._p['mg_noise'])

        rvec, vvec = gps.pull_gps(self.state, r_noise=self._p['gps_noise_r'], v_noise=self._p['gps_noise_v'])
        
        self.m.measurements.update(angularRate=wvec, quaternion=qvec, 
                                    bField=bvec, rvec=rvec, vvec=vvec)
        
        if self.mode == 'init':
            self.setpoint = np.array([-1, -1, -1, -1])
            self.q_dot = 0
        
        elif self.mode == 'detumble':
            self.setpoint = np.array([-1, -1, -1, -1])
            self.q_dot = 0

        elif self.mode == 'parking_orbit':
            self.setpoint = ground_target_guidance(self.m.measurements, 'nadir', [], 0)
            self.q_dot = 0
        
        elif 'transfer' in self.mode:
            gvec, q_dot = solveQLaw(self.state.orbit, self._p['Wp'], self._p['rp_min'], self._p['f_mag'],
                                    self._p['Wa'], self._p['We'], self._p['Wi'], self._p['Wraan'], self._p['Wargp'],
                                    self._p['aT'], self._p['eT'], self._p['iT'], self._p['raanT'], self._p['argpT'])
        
            r_mag = qm.normalize(self.m.measurements['rvec'])
            v_mag = qm.normalize(self.m.measurements['vvec'])
            h_mag = np.cross(v_mag, r_mag)
            
            gvec_local = [r_mag[0]*gvec[0] + v_mag[0]*gvec[1] + h_mag[0]*gvec[2],
                          r_mag[1]*gvec[0] + v_mag[1]*gvec[1] + h_mag[1]*gvec[2],
                          r_mag[2]*gvec[0] + v_mag[2]*gvec[1] + h_mag[2]*gvec[2]]
        
            self.setpoint = qm.axis_to_quaternion(np.array(gvec_local), 0)
            self.q_dot = q_dot
        
        elif 'hold' in self.mode:
            self.setpoint = ground_target_guidance(self.m.measurements, 'nadir', [], 0)
            self.q_dot = 0

    # control app
    # based on physical state of spacecraft and the mode, send commands
    # to the attitue actuators and thruster
    def _run_control(self):
        
        '''
        elif self.mode == 'transfer_1':
            if (self.state.orbit.nu << u.deg).value > 0 and (self.state.orbit.nu << u.deg).value < 180:
                if self.err < 0:
                    self._write_command_to_thruster(self._p['f_mag']*self.setpoint)
        '''

        thruster_command = 0
        rw_command = [0, 0, 0]
        mt_command = [0, 0, 0]
        q_error = [0, 0, 0, 1]
        

        if self.mode == 'detumble':
            mt_command = bDotController(self.m.measurements['angularRate'], self.m.measurements['bField'], self._p['bDot_gain'])


        elif self.mode == 'parking_orbit' or 'hold' in self.mode:
            q_error, rw_command = slidingModeController(self.m.measurements['angularRate'], 
                                                        self.m.measurements['quaternion'], self.setpoint, 
                                                        self._p['smc_kp'], self._p['smc_kd'], self._p['smc_sigma'], self._p['smc_order'])
        

        elif 'transfer' in self.mode:
            q_error, rw_command = slidingModeController(self.m.measurements['angularRate'], 
                                                        self.m.measurements['quaternion'], self.setpoint, 
                                                        self._p['smc_kp'], self._p['smc_kd'], self._p['smc_sigma'], self._p['smc_order'])
            
            if qm.magnitude([ q_error[0], q_error[1], q_error[2] ]) < 0.001:
                thruster_command = self._p['f_mag']

        
        self._write_command_to_thruster(thruster_command)
        self._write_command_to_magnetorquer(mt_command)
        self._write_command_to_reactionWheel(rw_command)

        self._write_variables_to_state(None, q_error)