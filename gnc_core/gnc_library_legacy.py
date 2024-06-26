import numpy as np
from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Earth
from astropy import units as u
import pandas as pd
import pyproj
from gnc_core.orbit_propagater import fpropagate
import matplotlib.pyplot as plt
import copy
import pyIGRF


from gnc_core.hardware_models.reaction_wheel import reactionWheelAssembly
from gnc_core.hardware_models.magnetorquers import magnetorquerAssembly


class satelliteState():
    def __init__(self, satellite_orbit, moment_of_inertia, satellite_mass,
                 satellite_area, q_initial, w_initial):
        self.I = moment_of_inertia
        self.m = satellite_mass
        self.A = satellite_area
        
        self.a = satellite_orbit[0] << u.m
        self.e = satellite_orbit[1] << u.one
        self.i = satellite_orbit[2] << u.deg
        self.O = satellite_orbit[3] << u.deg
        self.w = satellite_orbit[4] << u.deg
        self.nu = satellite_orbit[5] << u.deg
        
        self.orbit = Orbit.from_classical(Earth, self.a, self.e, self.i, self.O, self.w, self.nu)
        
        self.quaternion = q_initial
        self.angularRate = w_initial
        
        self.controlTorque = [0, 0, 0]
        self.disturbanceTorque = [0, 0, 0]
        self.controlForce = [0, 0, 0]
        self.disturbanceForce = [0, 0, 0]
        
        self.latitude = 0
        self.longitude = 0
        self.altitude = 0
        
        self.t = 0
        
    def propagateAttitude(self, controlTorque, disturbanceTorque, dt):
        
        self.controlTorque = controlTorque
        self.disturbanceTorque = disturbanceTorque
        
        w1 = copy.deepcopy(self.angularRate[0])
        w2 = copy.deepcopy(self.angularRate[1])
        w3 = copy.deepcopy(self.angularRate[2])
        
        q1 = copy.deepcopy(self.quaternion[0])
        q2 = copy.deepcopy(self.quaternion[1])
        q3 = copy.deepcopy(self.quaternion[2])
        q4 = copy.deepcopy(self.quaternion[3])
        
        omega_dot_1 = (1/self.I[0]) * ( (self.I[2] - self.I[1])*w2*w3 - (controlTorque[0] + disturbanceTorque[0]) )
        omega_dot_2 = (1/self.I[1]) * ( (self.I[0] - self.I[2])*w3*w1 - (controlTorque[1] + disturbanceTorque[1]) )
        omega_dot_3 = (1/self.I[2]) * ( (self.I[1] - self.I[0])*w1*w2 - (controlTorque[2] + disturbanceTorque[2]) )
        
        w1_new = w1 + omega_dot_1*dt
        w2_new = w2 + omega_dot_2*dt
        w3_new = w3 + omega_dot_3*dt
        
        w1_avg = 0.5 * (w1 + w1_new)
        w2_avg = 0.5 * (w2 + w2_new)
        w3_avg = 0.5 * (w3 + w3_new)
        
        q = quaternionIntegral([q1, q2, q3, q4], [w1_avg, w2_avg, w3_avg], dt)
        
        self.angularRate[0] = w1_new
        self.angularRate[1] = w2_new
        self.angularRate[2] = w3_new
        
        self.quaternion[0] = q[0]
        self.quaternion[1] = q[1]
        self.quaternion[2] = q[2]
        self.quaternion[3] = q[3]


        B_true = pyIGRF.igrf_value(self.latitude, self.longitude, self.altitude/10**3, 2024)
        B_true = (10**-9)*np.array( [B_true[3], B_true[4], B_true[5], 0] )


        B_body = quaternionMultiply( quaternionMultiply(q, B_true), conjugate(q) )
        self.B_true = np.array([ B_true[0], B_true[1], B_true[2] ])
        self.B_body = np.array([ B_body[0], B_body[1], B_body[2] ])

    
    def propagatePosition(self,  controlForce, disturbanceForce, dt):
        self.t += dt
        
        self.controlForce = controlForce
        self.disturbanceForce = disturbanceForce
        
        # new_orbit = self.orbit.propagate(dt << u.second)
        new_orbit = fpropagate(self.orbit, controlForce, self.m, self.A, self.t, dt)
        
        self.longitude, self.latitude, self.altitude = ECI_to_ECEF((new_orbit.r << u.meter).value, self.t)
        
        self.orbit = new_orbit
        
        


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
        quaternion_error = quaternionDifference(quaternion, setpoint)
        
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
        self.quaternion_error = quaternionDifference(self.quaternion, setpoint)
        
        controlTorque = []
        for i in range(3):
            if abs(self.quaternion_error[i]) >= self.sigma:
                u = self.kp * (self.quaternion_error[i] ** self.order) * sign(self.quaternion_error[3]) + self.kd * self.angularRate[i]
            
            else:
                u = self.kp * self.quaternion_error[i] * sign(self.quaternion_error[3]) + self.kd * self.angularRate[i]
                
            controlTorque.append(u)
            
        self.controlTorque = np.array( controlTorque )
        
        return self.controlTorque



def magnitude(vec):
    mag = 0
    for v in vec:
        mag += v**2
    
    return np.sqrt(mag)


def normalize(vec):
    
    m = magnitude(vec)
    
    v_ = []
    for v in vec:
        v_.append(v/m)
    
    return v_


def conjugate(q):
    return [-q[0], -q[1], -q[2], q[3]]



def quaternionIntegral(quaternion, angularRate, dt):
    
    wx, wy, wz = angularRate[0], angularRate[1], angularRate[2]
    qx, qy, qz, qw = quaternion[0], quaternion[1], quaternion[2], quaternion[3]
    
    return normalize( [ (0.5*dt) * ( wx*qw - wy*qz + wz*qy + qx ),
                        (0.5*dt) * ( wx*qz + wy*qw - wz*qx + qy ),
                        (0.5*dt) * (-wx*qy + wy*qx + wz*qw + qz ),
                        (0.5*dt) * (-wx*qx - wy*qy - wz*qz + qw ) ] )



def quaternionMultiply(q1, q2):
 
    w1, x1, y1, z1 = q1[3], q1[0], q1[1], q1[2]
    w2, x2, y2, z2 = q2[3], q2[0], q2[1], q2[2]
    
    return normalize( [ w1*x2 + w2*x1 + y1*z2 - y2*z1,
                        w1*y2 + w2*y1 - x1*z2 + x2*z1,
                        w1*z2 + w2*z1 + x1*y2 - x2*y1,
                        w1*w2 - x1*x2 - y1*y2 - z1*z2] )


def quaternionDifference(q1, q2):
    # https://math.stackexchange.com/questions/1782243/how-to-calculate-rotation-quaternion-between-two-orientation-quaternions
    # transforms from q1 to q2
    w1, x1, y1, z1 = q1[3], q1[0], q1[1], q1[2]
    w2, x2, y2, z2 = q2[3], q2[0], q2[1], q2[2]

    delta_w = w1*w2 - x1*x2 - y1*y2 - z1*z2
    delta_x = w1*x2 + w2*x1 - y1*z2 + y2*z1
    delta_y = w1*y2 + w2*y1 + x1*z2 - x2*z1
    delta_z = w1*z2 + w2*z1 - x1*y2 + x2*y1
    
    return normalize( [ delta_x, delta_y, delta_z, delta_w ] )


def quaternionDifferenceToAngularVelocity(q1, q2, dt):
    return (2 / dt) * np.array([ q1[0]*q2[1] - q1[1]*q2[0] - q1[2]*q2[3] + q1[3]*q2[2],
                                 q1[0]*q2[2] + q1[1]*q2[3] - q1[2]*q2[0] - q1[3]*q2[1],
                                 q1[0]*q2[3] - q1[1]*q2[2] + q1[2]*q2[1] - q1[3]*q2[0] ])



def ECI_to_ECEF(r, t):
    
    gamma = (360/(23.9345*3600)) * t
    gamma *= np.pi/180
    
    x = r[0]*np.cos(gamma) - r[1]*np.sin(gamma) # meters
    y = r[0]*np.sin(gamma) + r[1]*np.cos(gamma) # meters
    z = r[2] # meters
    
    transformer = pyproj.Transformer.from_crs({"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
                                              {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'})
    
    lon, lat, alt = transformer.transform(x, y, z, radians=False)
    
    return lon, lat, alt


def axis_to_quaternion(axis, rotation):
    return normalize( np.array([ axis[0]*np.sin(rotation/2),
                                 axis[1]*np.sin(rotation/2),
                                 axis[2]*np.sin(rotation/2),
                                 np.cos(rotation/2) ]) )


def sign(val):
    if val == 0: return 0
    elif val < 0: return -1
    else: return 1
    

def appendDataFrame(df, state, t, error, rwheel):
    tempdf = pd.DataFrame({
                           'time':[t], 
                           'Mx':[state.controlTorque[0]], 'My':[state.controlTorque[1]], 'Mz':[state.controlTorque[2]],
                           'Fx':[state.controlForce[0]], 'Fy':[state.controlForce[1]], 'Fz':[state.controlForce[2]],
                           'wx':[state.angularRate[0]], 'wy':[state.angularRate[1]], 'wz':[state.angularRate[2]],
                           'q1':[state.quaternion[0]], 'q2':[state.quaternion[1]], 'q3':[state.quaternion[2]], 'q4':[state.quaternion[3]],
                           'q1e':[error[0]], 'q2e':[error[1]], 'q3e':[error[2]], 'q4e':[error[3]],
                           'wheel1':[rwheel.wheel_speeds[0]], 'wheel2':[rwheel.wheel_speeds[1]], 'wheel3':[rwheel.wheel_speeds[2]],
                           'a':[state.orbit.a.value], 'e':[state.orbit.ecc.value], 'i':[state.orbit.inc.value], 
                           'raan':[state.orbit.raan.value], 'argp':[state.orbit.argp.value], 'nu':[state.orbit.nu.value],
                           'x':[(state.orbit.r[0] << u.km).value], 'y':[(state.orbit.r[1] << u.km).value], 'z':[(state.orbit.r[2] << u.km).value],
                           'vx':[(state.orbit.v[0] << u.km/u.second).value], 'vy':[(state.orbit.v[1] << u.km/u.second).value], 'vz':[(state.orbit.v[2] << u.km/u.second).value],
                           'Bbx':[state.B_body[0]], 'Bby':[state.B_body[1]], 'Bbz':[state.B_body[2]],
                           'Btx':[state.B_true[0]], 'Bty':[state.B_true[1]], 'Btz':[state.B_true[2]],
                           'latitude':[state.latitude], 'longitude':[state.longitude], 'altitude':[state.altitude]
                           })
    return pd.concat([df, tempdf])



def plot_ground_track(df):
    img = plt.imread("Earth.jpg")
    fig, ax = plt.subplots()
    ax.imshow(img, extent=[-180, 180, -90, 90])
    plt.scatter(np.array(df[df['longitude'] >= 0]['longitude']), np.array(df[df['longitude'] >= 0]['latitude']), color='yellow', s=0.5)
    plt.scatter(np.array(df[df['longitude'] <= 0]['longitude']), np.array(df[df['longitude'] <= 0]['latitude']), color='yellow', s=0.5)

    plt.show()
    

def plot_orbit(df):
    
    x = np.array(df['x'])
    y = np.array(df['y'])
    z = np.array(df['z'])
    
    ax = plt.figure().add_subplot(projection='3d')
    ax.plot(x, y, z, lw=0.5)
    
    for i in range(0, len(df), int(1 + len(df)/50)):
        ux = 100*df['Fx'].iloc[i]
        uy = 100*df['Fy'].iloc[i]
        uz = 100*df['Fz'].iloc[i]
        ax.quiver(x[i], y[i], z[i],   ux, uy, uz)
    
    v = []
    for i in range(len(df)):
        v.append( magnitude([x[i], y[i], z[i]]) )
    
    limit = max(v)
    
    ax.set_xlim([-limit-5000, limit+5000])
    ax.set_ylim([-limit-5000, limit+5000])
    ax.set_zlim([-limit-5000, limit+5000])

    plt.show()


def plot_quaternion_error(df):
    t = np.array( df['time'] )
    
    q1e = np.array( df['q1e'] )
    q2e = np.array( df['q2e'] )
    q3e = np.array( df['q3e'] )
    q4e = np.array( df['q4e'] )
    
    q1 = np.array( df['q1'] )
    q2 = np.array( df['q2'] )
    q3 = np.array( df['q3'] )
    q4 = np.array( df['q4'] )

    ################

    plt.figure(1)
    plt.plot(t, q1e, label='q1e')
    plt.plot(t, q2e, label='q2e')
    plt.plot(t, q3e, label='q3e')
    plt.plot(t, q4e, label='q4e')
    
    plt.grid()
    plt.legend()
    
    plt.xlabel('time (s)')
    plt.ylabel('quaternion error')

    ################

    plt.figure(2)
    plt.plot(t, q1, label='q1')
    plt.plot(t, q2, label='q2')
    plt.plot(t, q3, label='q3')
    plt.plot(t, q4, label='q4')
    
    plt.grid()
    plt.legend()
    
    plt.xlabel('time (s)')
    plt.ylabel('quaternion')

    ################

    plt.figure(3)
    plt.plot(t, df['wheel1'], label='wheel 1 speed')
    plt.plot(t, df['wheel2'], label='wheel 2 speed')
    plt.plot(t, df['wheel3'], label='wheel 3 speed')
    
    plt.grid()
    plt.legend()
    
    plt.xlabel('time (s)')
    plt.ylabel('wheel speeds (rad/s)')

    ################

    plt.figure(4)
    plt.plot(t, df['wx'], label='wx')
    plt.plot(t, df['wy'], label='wy')
    plt.plot(t, df['wz'], label='wz')
    
    plt.grid()
    plt.legend()
    
    plt.xlabel('time (s)')
    plt.ylabel('angular rate (rad/s)')

    ################

    plt.show()

