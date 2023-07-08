import numpy as np
from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Earth
from astropy import units as u
import pandas as pd
import pyproj
from orbit_propagater import fpropagate
import matplotlib.pyplot as plt



class reactionWheelAssembly():
    def __init__(self, wheel_moi, max_wheel_speed, max_wheel_torque, min_wheel_torque):
        self.wheel_moi = wheel_moi
        self.max_wheel_speed = max_wheel_speed
        self.max_wheel_torque = max_wheel_torque
        
        self.wheel_speeds = np.array([0, 0, 0])
        self.wheel_torques = np.array([0, 0, 0])
        
        self.torque_command = [0, 0, 0]
    
    def actuateReactionWheels(self, torque_command, timestep):
        for i in range(3):
            self.torque_command[i] = limit(torque_command[i], self.max_wheel_torque)
        
        wheel_speeds_next = []
        wheel_torques_next = []
        for i in range(3):
            w = self.wheel_speeds[i] + (1 / self.wheel_moi) * self.torque_command[i] * timestep

            if abs(w) < self.max_wheel_speed:
                wheel_torques_next.append( self.torque_command[i] )
                wheel_speeds_next.append( w )
            else:
                wheel_torques_next.append( 0 )
                wheel_speeds_next.append( self.wheel_speeds[i] )
        
        self.wheel_torques = wheel_torques_next
        self.wheel_speeds = wheel_speeds_next



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
        
        w1 = self.angularRate[0]
        w2 = self.angularRate[1]
        w3 = self.angularRate[2]
        
        q1 = self.quaternion[0]
        q2 = self.quaternion[1]
        q3 = self.quaternion[2]
        q4 = self.quaternion[3]
        
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
        
    
    def propagatePosition(self,  controlForce, disturbanceForce, dt):
        self.t += dt
        
        self.controlForce = controlForce
        self.disturbanceForce = disturbanceForce
        
        new_orbit = self.orbit.propagate(dt << u.second)
        new_orbit = fpropagate(new_orbit, controlForce, self.m, self.A, self.t, dt)
        
        self.longitude, self.latitude, self.altitude = ECI_to_ECEF((new_orbit.r << u.meter).value, self.t)
        
        self.orbit = new_orbit
        
        


class p2Controller():
    def __init__(self, kp, kd):
        self.kp = kp
        self.kd = kd
    
    def input_function(self, angularRate, quaternion):
        self.angularRate = angularRate
        self.quaternion = quaternion
    
    def output_function(self, setpoint):
        self.quaternion_error = quaternionMultiply( self.quaternion, conjugate(setpoint) )
        
        controlTorque = []
        for i in range(3):
            u = self.kp * self.quaternion_error[i] * self.quaternion_error[3] + self.kd * self.angularRate[i]
            controlTorque.append(u)
        
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
        self.quaternion_error = quaternionMultiply( self.quaternion, conjugate(setpoint) )
        
        controlTorque = []
        for i in range(3):
            if abs(self.quaternion_error[i]) >= self.sigma:
                u = self.kp * (self.quaternion_error[i] ** self.order) * self.quaternion_error[3] + self.kd * self.angularRate[i]
            
            else:
                u = self.kp * self.quaternion_error[i] * self.quaternion_error[3] + self.kd * self.angularRate[i]
                
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
    n = np.array([0, 0, 0])
    
    if type(q1) != type(n):
        q1 = np.array(q1)
        
    if type(q2) != type(n):
        q2 = np.array(q2)
    
    
    w1, x1, y1, z1 = q1[3], q1[0], q1[1], q1[2]
    w2, x2, y2, z2 = q2[3], q2[0], q2[1], q2[2]
    
    return normalize( [ w1*x2 + w2*x1 + y1*z2 - y2*z1,
                        w1*y2 + w2*y1 - x1*z2 + x2*z1,
                        w1*z2 + w2*z1 + x1*y2 - x2*y1,
                        w1*w2 - x1*x2 - y1*y2 - z1*z2] )


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





def limit(val, maximum):
    if val > maximum:       return maximum
    elif val < -maximum:    return -maximum
    else:                   return val


def sign(val):
    if val == 0: return 0
    elif val < 0: return -1
    else: return 1
    

def appendDataFrame(df, state, t, error):
    tempdf = pd.DataFrame({
                           'time':[t], 
                           'Mx':[state.controlTorque[0]], 'My':[state.controlTorque[1]], 'Mz':[state.controlTorque[2]],
                           'Fx':[state.controlForce[0]], 'Fy':[state.controlForce[1]], 'Fz':[state.controlForce[2]],
                           'w1':[state.angularRate[0]], 'w2':[state.angularRate[1]], 'w3':[state.angularRate[2]],
                           'q1':[state.quaternion[0]], 'q2':[state.quaternion[1]], 'q3':[state.quaternion[2]], 'q4':[state.quaternion[3]],
                           'q1e':[error[0]], 'q2e':[error[1]], 'q3e':[error[2]], 'q4e':[error[3]],
                           'a':[state.orbit.a.value], 'e':[state.orbit.ecc.value], 'i':[state.orbit.inc.value], 
                           'raan':[state.orbit.raan.value], 'argp':[state.orbit.argp.value], 'nu':[state.orbit.nu.value],
                           'x':[(state.orbit.r[0] << u.km).value], 'y':[(state.orbit.r[1] << u.km).value], 'z':[(state.orbit.r[2] << u.km).value],
                           'latitude':[state.latitude], 'longitude':[state.longitude], 'altitude':[state.altitude]
                           })
    return pd.concat([df, tempdf])



def plot_ground_track(df):
    img = plt.imread("Earth.jpg")
    fig, ax = plt.subplots()
    ax.imshow(img, extent=[-180, 180, -90, 90])
    plt.scatter(np.array(df[df['longitude'] >= 0]['longitude']), np.array(df[df['longitude'] >= 0]['latitude']), color='yellow', s=0.5)
    plt.scatter(np.array(df[df['longitude'] <= 0]['longitude']), np.array(df[df['longitude'] <= 0]['latitude']), color='yellow', s=0.5)

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

def plot_quaternion_error(df):
    t = np.array( df['time'] )
    
    q1e = np.array( df['q1e'] )
    q2e = np.array( df['q2e'] )
    q3e = np.array( df['q3e'] )
    q4e = np.array( df['q4e'] )
    
    plt.figure(1)
    plt.plot(t, q1e, label='q1e')
    plt.plot(t, q2e, label='q2e')
    plt.plot(t, q3e, label='q3e')
    plt.plot(t, q4e, label='q4e')
    
    plt.grid()
    plt.legend()
    
    plt.xlabel('time (s)')
    plt.ylabel('quaternion error')