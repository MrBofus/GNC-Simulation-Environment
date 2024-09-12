import quaternion_math.quaternionMath as qm

import gnc_core.gnc_library as gnc
import numpy as np



def sign(val):
    if val == 0: return 0
    elif val < 0: return -1
    else: return 1


def p2controller(angularRate, quaternion, setpoint, kp, kd):
    # https://ltu.diva-portal.org/smash/get/diva2:1010947/FULLTEXT01.
    
    quaternion_error = qm.quaternionDifference(quaternion, setpoint)
        
    controlTorque = []
    for i in range(3):
        u = ( kp * quaternion_error[i] * sign(quaternion_error[3]) ) + ( kd * angularRate[i] )
        controlTorque.append(u)

    return quaternion_error, np.array( controlTorque )


def slidingModeController(angularRate, quaternion, setpoint,
                            kp, kd, sigma, order):
    if setpoint == [-1, -1, -1, -1]:
        return [0, 0, 0, 1], np.array([0, 0, 0])
    
    quaternion_error = qm.quaternionDifference(quaternion, setpoint)
        
    controlTorque = []
    for i in range(3):
        if abs(quaternion_error[i]) >= sigma:
            u = kp * (quaternion_error[i] ** order) * sign(quaternion_error[3]) + kd * angularRate[i]
        
        else:
            u = kp * quaternion_error[i] * sign(quaternion_error[3]) + kd * angularRate[i]
            
        controlTorque.append(u)

    return quaternion_error, np.array( controlTorque )


'''
def magneticController(angularRate, quaternion, setpoint, quaternion_error_previous,
                            kp, ki, kd, sigma, order, magnetic_field):
    if setpoint == [-1, -1, -1, -1]:
        return [0, 0, 0, 1], np.array([0, 0, 0])

    quaternion_error = qm.quaternionDifference(quaternion, setpoint)

    controlTorque = []
    for i in range(3):

        if abs(quaternion_error[i]) >= sigma:
            u = kp[i] * (quaternion_error[i] ** order) * sign(quaternion_error[3]) + kd[i] * angularRate[i]
        
        else:
            u = kp[i] * quaternion_error[i] * sign(quaternion_error[3]) + kd[i] * angularRate[i]
        
        # print(kp[i] * quaternion_error[i] * sign(quaternion_error[3]), kd[i] * angularRate[i])
        controlTorque.append(u)
    # print('--')
    controlTorque = -np.cross(controlTorque, qm.normalize(magnetic_field))
    return quaternion_error, controlTorque
'''

def magneticController(angularRate, quaternion, 
                       q_setpoint, w_setpoint, q_err_prev,
                       kp, ki, kd, sigma, order, magnetic_field):
    
    if q_setpoint == [-1, -1, -1, -1]:
        return [0, 0, 0, 1], np.array([0, 0, 0])

    quaternion_error = qm.quaternionDifference(quaternion, q_setpoint)

    controlTorque = []
    for i in range(3):


        if abs(quaternion_error[i]) >= sigma:
            u = kp[i] * (quaternion_error[i] ** order) + ki[i] * q_err_prev[i] + kd[i] * (angularRate[i] - w_setpoint[i])
        
        else:
            u = kp[i] * quaternion_error[i] + ki[i] * q_err_prev[i] + kd[i] * (angularRate[i] - w_setpoint[i])
 
        controlTorque.append(u)


    controlTorque = -np.cross(controlTorque, qm.normalize(magnetic_field))

    return quaternion_error, controlTorque


def bDotController(angularRate, bField, gain):
    return -gain*np.cross(angularRate, qm.normalize(bField))


def ground_target_guidance(measurements, mode, ground_station_list, gs_range):
    rvec = np.array( qm.normalize(measurements['rvec']) )
    vvec = np.array( qm.normalize(measurements['vvec']) )
    hvec = np.array( np.cross(rvec, vvec) )
    
    if mode == 'nadir':
        return qm.dcm_to_quaternion( np.array([hvec, vvec, -rvec]) )
    
    elif mode == 'prograde':
        return qm.dcm_to_quaternion( np.array([hvec, rvec, vvec]) )
    
    elif mode == 'downlink':
        for i in range(len(ground_station_list)):
            g_ = gnc.ECEF_to_ECI( ground_station_list[i], 0 )

            r_ = 1000*np.array(measurements['rvec']) - np.array(g_)
            if qm.magnitude(r_) < gs_range:
                v_ = qm.orthogonalize(measurements['vvec'], r_)
                h_ = np.cross(r_, v_)

                return qm.dcm_to_quaternion(np.array(np.array([h_, v_, -r_])))  
        return qm.dcm_to_quaternion( np.array([hvec, vvec, -rvec]) )


def checkVicinity(measurements, ground_station_list, gs_range):
    for i in range(len(ground_station_list)):
        g_ = gnc.ECEF_to_ECI( ground_station_list[i], 0 )
        r_ = 1000*np.array(measurements['rvec']) - np.array(g_)

        if qm.magnitude(r_) < gs_range:
            return r_

    return [-1, -1, -1]




def printout(str):
    print('\033[0;32m  from flight_software.schedulerApp._update_user_variables():\n' + str + '\033[0m')
