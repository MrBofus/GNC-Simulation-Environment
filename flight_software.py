import quaternionMath as qm
import numpy as np
import copy

import gnc_core.hardware_models.magnetometer as mm
import gnc_core.hardware_models.rate_gyro as rg
import gnc_core.hardware_models.star_tracker as st
import gnc_core.hardware_models.gps as gps



def sign(val):
    if val == 0: return 0
    elif val < 0: return -1
    else: return 1

'''
class p2Controller():
    # https://ltu.diva-portal.org/smash/get/diva2:1010947/FULLTEXT01.pdf
    def __init__(self, kp, kd):
        self.kp = kp
        self.kd = kd
    
    def input_function(self, angularRate, quaternion):
        self.angularRate = angularRate
        self.quaternion = quaternion
    
    def output_function(self, setpoint):
        quaternion, angularRate = copy.deepcopy(self.quaternion), copy.deepcopy(self.angularRate)
        # quaternion_error = quaternionMultiply( quaternion, conjugate(setpoint) )
        quaternion_error = qm.quaternionDifference(quaternion, setpoint)
        
        controlTorque = []
        for i in range(3):
            u = ( self.kp * quaternion_error[i] * sign(quaternion_error[3]) ) + ( self.kd * angularRate[i] )
            controlTorque.append(u)

        self.quaternion_error = quaternion_error
        self.controlTorque = np.array( controlTorque )
        
        return self.controlTorque



class smcController():
    def __init__(self, kp, kd, sigma, order):
        self.kp = kp
        self.kd = kd
        
        self.sigma = sigma
        self.order = order
    
    def input_function(self, angularRate, quaternion):
        self.angularRate = angularRate
        self.quaternion = quaternion
    
    def output_function(self, setpoint):
        self.quaternion_error = qm.quaternionDifference(self.quaternion, setpoint)
        
        controlTorque = []
        for i in range(3):
            if abs(self.quaternion_error[i]) >= self.sigma:
                u = self.kp * (self.quaternion_error[i] ** self.order) * sign(self.quaternion_error[3]) + self.kd * self.angularRate[i]
            
            else:
                u = self.kp * self.quaternion_error[i] * sign(self.quaternion_error[3]) + self.kd * self.angularRate[i]
                
            controlTorque.append(u)
            
        self.controlTorque = np.array( controlTorque )
        
        return self.controlTorque
'''

def p2controller(angularRate, quaternion, setpoint, kp, kd):
    quaternion_error = qm.quaternionDifference(quaternion, setpoint)
        
    controlTorque = []
    for i in range(3):
        u = ( kp * quaternion_error[i] * sign(quaternion_error[3]) ) + ( kd * angularRate[i] )
        controlTorque.append(u)

    return quaternion_error, np.array( controlTorque )


def slidingModeController(angularRate, quaternion, setpoint,
                            kp, kd, sigma, order):
    quaternion_error = qm.quaternionDifference(quaternion, setpoint)
        
    controlTorque = []
    for i in range(3):
        if abs(quaternion_error[i]) >= sigma:
            u = kp * (quaternion_error[i] ** order) * sign(quaternion_error[3]) + kd * angularRate[i]
        
        else:
            u = kp * quaternion_error[i] * sign(quaternion_error[3]) + kd * angularRate[i]
            
        controlTorque.append(u)

    return quaternion_error, np.array( controlTorque )


def bDotController(angularRate, bField, gain):
    return -gain*np.cross(angularRate, bField)



def printout(str):
    print('\033[0;32m\tfrom flight_software.schedulerApp():\n' + str + '\033[0m')



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
                    'bDot_gain': 0.1,
            
                    'smc_kp':1.0,
                    'smc_kd':10.0,
                    'smc_sigma':0.1,
                    'smc_order':3,

                    'rg_noise':1/5000,
                    'rg_bias':0.001,

                    'st_noise':1/5000,
                    
                    'mg_noise':10**-11,
                    
                    'gps_noise_r':0.01,
                    'gps_noise_v':0.01,

                    }
        
        self._p = defaults
        outstr = ''
        for arg, value in kwargs.items():
            if arg in defaults:
                self._p[arg]=value
                outstr += '\t--- updated: ' + arg + ' --> ' + str(value) + '\n'
        
        outstr += '\n' + str(self._p)
        printout(outstr)


    # interface functions {#121, 14}
    def _physics_to_harware_int(self, state):
        self._s = state
        self.systemClock += 1

    def _write_command_to_reactionWheel(self, command):
        self.rw_command = command
    
    def _write_command_to_magnetorquer(self, command):
        self.mt_command = command
    
    def _write_variables_to_state(self, mode, q_error):
        self.mode = mode
        self.q_error = q_error
    

    # main function {#112, 9}
    def _iterate(self):
        if self.systemClock == 1:
            self.mode = 'detumble'
        
        self._run_navigation()
        self._run_control()

        pass


    # flight apps {#211, 45}
    def _run_navigation(self):
        wvec = rg.pull_gyro(self._s, noise=self._p['rg_noise'])
        qvec = st.pull_star_tracker(self._s, noise=self._p['st_noise'])
        bvec = mm.pull_magnetometer(self._s, noise=self._p['mg_noise'])

        rvec, vvec = gps.pull_gps(self._s, r_noise=self._p['gps_noise_r'], v_noise=self._p['gps_noise_v'])
        
        self.m.measurements.update(angularRate=wvec, quaternion=qvec, 
                                    bField=bvec, rvec=rvec, vvec=vvec)

    def _run_control(self):
        if self.mode == 'detumble':
            mt_command = bDotController(self.m.measurements['angularRate'], self.m.measurements['bField'], self._p['bDot_gain'])
            rw_command = [0, 0, 0]
            q_error = [0, 0, 0, 1]

            if qm.magnitude(self.m.measurements['angularRate']) < 0.002:
                self._write_command_to_magnetorquer(mt_command)
                self._write_command_to_reactionWheel(rw_command)
                self._write_variables_to_state('dataCollection', q_error)
                return None
            
            self._write_command_to_magnetorquer(mt_command)
            self._write_command_to_reactionWheel(rw_command)
            self._write_variables_to_state(self.mode, q_error)
            return None
        
        elif self.mode == 'dataCollection':
            mt_command = [0, 0, 0]

            rvec = np.array( qm.normalize(self.m.measurements['rvec']) )
            vvec = np.array( qm.normalize(self.m.measurements['vvec']) )
            hvec = np.array( np.cross(rvec, vvec) )

            # nadir_direction = qm.dcm_to_quaternion( np.array([hvec, -vvec, rvec]) )
            nadir_direction = qm.dcm_to_quaternion( np.array([hvec, vvec, -rvec]) )

            q_error, rw_command = slidingModeController(self.m.measurements['angularRate'], 
                                                        self.m.measurements['quaternion'], nadir_direction, 
                                                        self._p['smc_kp'], self._p['smc_kd'], self._p['smc_sigma'], self._p['smc_order'])
            
            self._write_command_to_magnetorquer(mt_command)
            self._write_command_to_reactionWheel(rw_command)
            self._write_variables_to_state(self.mode, q_error)
            return None