import quaternion_math.quaternionMath as qm

import gnc_core.gnc_library as gnc
import numpy as np



class slidingModeController():
    def __init__(self, kp, kd, window):
        self.kp = kp
        self.kd = kd
        self.window = window
        self.order = 3

    def run(self, q_setpoint, w_setpoint, measurements):
        if q_setpoint[0] == -1 and q_setpoint[1] == -1 and q_setpoint[2] == -1:
            return [0, 0, 0, 1], np.array([0, 0, 0])
        
        q = np.array(measurements['quaternion'])
        w = np.array(measurements['angularRate'])

        q_error = qm.quaternionDifference(q, q_setpoint)
        w_error = w - w_setpoint

        # s = (w - w_setpoint) + self.lambd * q_error_vec * np.sign(q_error[3])
        
        u = []
        for i in range(3):
            if abs(q_error[i]) >= self.window:
                u_ = self.kp * (q_error[i] ** 3) * np.sign(q_error[3]) + self.kd * w_error[i]
            else:
                u_ = self.kp * q_error[i] * np.sign(q_error[3]) + self.kd * w_error[i]
            u.append(u_)

        return q_error, np.array(u)


def bDotController(angularRate, bField, gain):
    return -gain*np.cross(angularRate, bField)


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
