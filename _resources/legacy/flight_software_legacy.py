import quaternionMath as qm
import numpy as np
import copy

import gnc_core.hardware_models.magnetometer as mm
import gnc_core.hardware_models.rate_gyro as rg
import gnc_core.hardware_models.star_tracker as st
import gnc_core.hardware_models.gps as gps


kp = 1.0
kd = 10.0
sigma = 0.1
order = 3


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

        

def _navigation(state, measurements):
    measured_angular_rate = rg.pull_gyro(state, noise=(1/5000))
    measured_quaternion = st.pull_star_tracker(state, noise=(1/5000))
    measured_Bfield = mm.pull_magnetometer(state, noise=10**-11)

    rvec, vvec = gps.pull_gps(state, r_noise=0.01, v_noise=0.01)
    
    measurements.update(measured_angular_rate, measured_quaternion, 
                        measured_Bfield, rvec, vvec)



def _control_detumble(measurements):
    mt_command = bDotController(measurements.angularRate, measurements.bField, 10**-1)
    rw_command = [0, 0, 0]
    q_error = [0, 0, 0, 1]

    if qm.magnitude(measurements.angularRate) < 0.001:
        return mt_command, rw_command, 'operation', q_error
    
    return mt_command, rw_command, 'detumble', q_error


def _control_dataCollection(measurements):
    mt_command = [0, 0, 0]

    rvec = qm.normalize(measurements.r)
    vvec = qm.normalize(measurements.v)

    hvec = np.cross(rvec, vvec)

    _ = qm.axis_to_quaternion(vvec, np.pi/2)

    q_error, rw_command = slidingModeController(measurements.angularRate, 
                                                measurements.quaternion, 
                                                [0, 0, 0, 1], kp, kd, sigma, order)
    
    return mt_command, rw_command, 'dataCollection', q_error



class schedulerApp():
    class measured_state():
        def __init__(self):
            self.angularRate = np.array([0, 0, 0])
            self.quaternion = np.array([0, 0, 0])
            self.bField = np.array([0, 0, 0])
            self.r = np.array([0, 0, 0])
            self.v = np.array([0, 0, 0])

        def update(self, measured_angular_rate, measured_quaternion, 
                    measured_Bfield, rvec, vvec):
            self.angularRate = measured_angular_rate
            self.quaternion = measured_quaternion
            self.bField = measured_Bfield
            self.r = rvec
            self.v = vvec
            
    def __init__(self):
        pass

    def _update_user_variables(self):
        pass

    def _run_navigation(self):
        pass

    def _run_control(self):
        pass

    def _update_global_state(self):
        pass