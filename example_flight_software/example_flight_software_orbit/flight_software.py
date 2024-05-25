from example_flight_software.example_flight_software_orbit.QLaw import solveQLaw

import copy
import quaternion_math.quaternionMath as qm
from astropy import units as u
import numpy as np


class schedulerApp():
    # measurements class {#4f7, 9}
    class measured_state():
        def __init__(self):
            self.measurements = {}

        def update(self, **kwargs):
            self.measurements = {}
            for arg, value in kwargs.items():
                self.measurements.update({arg:value})


    # init parameters {#121, 5}
    def __init__(self):
        self._p = {}
        self.m = self.measured_state()
        self.systemClock = 0


    # user-defined variables {#121, 31}
    def _update_user_variables(self, **kwargs):
        defaults = {

                    'Wp':1,
                    'rp_min':400000,

                    'f_mag':330 * 10**-6,

                    'Wa':10,
                    'We':10,
                    'Wi':1,
                    'Wraan':0,
                    'Wargp':0,

                    'aT':(530+6378.1)*10**3,
                    'eT':0.001,
                    'iT':51.9 * np.pi/180,
                    'raanT':0,
                    'argpT':0,

                    }
        
        self._p = defaults
        outstr = ''
        for arg, value in kwargs.items():
            if arg in defaults:
                outstr += '\t--- updated: ' + arg + '  ' + str(defaults[arg]) + ' --> ' + str(value) + '\n'
                self._p[arg]=value
            else:
                outstr += '\033[0;31m\t--- warning: ' + arg + ' not recognized as defualt variable, no action taken\033[0;32m\n'
        
        outstr += '\n\t' + str(self._p)
        print('\n\nscheduler:\n'+outstr)


    def _write_command_to_thruster(self, command):
        self.thruster_command = command

    def _write_variables_to_state(self, mode):
        self.mode = mode
    

    # main function {#112, 9}
    def _iterate(self, state):
        self.state = copy.deepcopy(state)
        self.systemClock += 1
        
        if self.systemClock == 1:
            self.mode = 'init'
        
        self._run_guidance()
        self._run_navigation()
        self._run_control()

        pass


    # flight apps {#211, 48}
    def _run_guidance(self):
        if self.mode == 'init':
            self.setpoint = np.array([0, 0, 0])
            pass
        
        else:
            gvec, _ = solveQLaw(self.state.orbit, self._p['Wp'], self._p['rp_min'], self._p['f_mag'],
                                    self._p['Wa'], self._p['We'], self._p['Wi'], self._p['Wraan'], self._p['Wargp'],
                                    self._p['aT'], self._p['eT'], self._p['iT'], self._p['raanT'], self._p['argpT'])
        
            r_mag = qm.normalize(self.state.orbit.r.value)
            v_mag = qm.normalize(self.state.orbit.v.value)
            h_mag = np.cross(v_mag, r_mag)
            
            gvec_local = [r_mag[0]*gvec[0] + v_mag[0]*gvec[1] + h_mag[0]*gvec[2],
                          r_mag[1]*gvec[0] + v_mag[1]*gvec[1] + h_mag[1]*gvec[2],
                          r_mag[2]*gvec[0] + v_mag[2]*gvec[1] + h_mag[2]*gvec[2]]
        
            self.setpoint = np.array(gvec_local)


    def _run_navigation(self):
        if self.mode == 'init':

            if self.systemClock > 2* (24*3600)/(3*60):
                self._write_variables_to_state('transfer_1')
                pass
            
            pass
        
        elif self.mode == 'transfer_1':

            if (self.state.orbit.a << u.meter).value > 0.9999*self._p['aT'] and (self.state.orbit.a << u.meter).value < 1.0001*self._p['aT']:
                if (self.state.orbit.inc << u.radian).value > 0.9999*self._p['iT'] and (self.state.orbit.inc << u.radian).value < 1.0001*self._p['iT']:
                    if self.state.orbit.ecc.value > 0.9999*self._p['eT'] and self.state.orbit.ecc.value < 1.0001*self._p['eT']:
                        self._write_variables_to_state('hold_1')
                        self.t_crit = self.systemClock

                        pass
            pass
        
        elif self.mode == 'hold_1':

            if self.systemClock > self.t_crit + 2* (24*3600)/(3*60):
                self._write_variables_to_state('transfer_2')
                self._update_user_variables(aT=self._p['aT']+200000, Wa=15, We=1, Wi=1)
                
                pass
            pass
        
        elif self.mode == 'transfer_2':

            if (self.state.orbit.a << u.meter).value > 0.9999*self._p['aT'] and (self.state.orbit.a << u.meter).value < 1.0001*self._p['aT']:
                if (self.state.orbit.inc << u.radian).value > 0.9999*self._p['iT'] and (self.state.orbit.inc << u.radian).value < 1.0001*self._p['iT']:
                    if self.state.orbit.ecc.value > 0.9999*self._p['eT'] and self.state.orbit.ecc.value < 1.01*self._p['eT']:
                        self._write_variables_to_state('hold_2')
                        self.t_crit = self.systemClock

                        pass
            pass
        
        elif self.mode == 'hold_2':

            if self.systemClock > self.t_crit + 2* (24*3600)/(3*60):
                self._write_variables_to_state('exit')
                pass
            pass


    def _run_control(self):
        
        if self.mode == 'init':
            self._write_command_to_thruster([0, 0, 0])
            pass
        
        elif self.mode == 'transfer_1':
            self._write_command_to_thruster(self._p['f_mag']*self.setpoint)
            pass
        
        elif self.mode == 'hold_1':
            self._write_command_to_thruster([0, 0, 0])
            pass
        
        elif self.mode == 'transfer_2':
            self._write_command_to_thruster(self._p['f_mag']*self.setpoint)
            pass
        
        elif self.mode == 'hold_2':
            self._write_command_to_thruster([0, 0, 0])
            pass