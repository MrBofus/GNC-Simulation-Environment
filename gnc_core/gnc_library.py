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
import quaternion_math.quaternionMath as qm

import warnings
warnings.filterwarnings('ignore')

from gnc_core.hardware_models.reaction_wheel import reactionWheelAssembly
from gnc_core.hardware_models.magnetorquers import magnetorquerAssembly

import numba as nb

@nb.jit(nopython=True)
def _fast_integrator(controlTorque, disturbanceTorque, angularRate, quaternion, I, dt):
    
    w1 = angularRate[0]
    w2 = angularRate[1]
    w3 = angularRate[2]
    
    q1 = quaternion[0]
    q2 = quaternion[1]
    q3 = quaternion[2]
    q4 = quaternion[3]
    
    omega_dot_1 = (1/I[0]) * ( (I[2] - I[1])*w2*w3 - (controlTorque[0] + disturbanceTorque[0]) )
    omega_dot_2 = (1/I[1]) * ( (I[0] - I[2])*w3*w1 - (controlTorque[1] + disturbanceTorque[1]) )
    omega_dot_3 = (1/I[2]) * ( (I[1] - I[0])*w1*w2 - (controlTorque[2] + disturbanceTorque[2]) )
    
    w1_new = w1 + omega_dot_1*dt
    w2_new = w2 + omega_dot_2*dt
    w3_new = w3 + omega_dot_3*dt
    
    w1_avg = 0.5 * (w1 + w1_new)
    w2_avg = 0.5 * (w2 + w2_new)
    w3_avg = 0.5 * (w3 + w3_new)
    
    q = qm.quaternionIntegral([q1, q2, q3, q4], [w1_avg, w2_avg, w3_avg], dt)
    
    angularRate[0] = w1_new
    angularRate[1] = w2_new
    angularRate[2] = w3_new
    
    quaternion[0] = q[0]
    quaternion[1] = q[1]
    quaternion[2] = q[2]
    quaternion[3] = q[3]
    
    return angularRate, quaternion, controlTorque, disturbanceTorque



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
        
        self.angularRate, self.quaternion, self.controlTorque, self.disturbanceTorque = _fast_integrator(controlTorque, disturbanceTorque, self.angularRate, self.quaternion, self.I, dt)


        B_true = pyIGRF.igrf_value(self.latitude, self.longitude, self.altitude/10**3, 2024)
        B_true = (10**-9)*np.array( [B_true[3], B_true[4], B_true[5], 0] )


        B_body = qm.quaternionMultiply( qm.quaternionMultiply(self.quaternion, B_true), qm.conjugate(self.quaternion) )
        self.B_true = np.array([ B_true[0], B_true[1], B_true[2] ])
        self.B_body = qm.magnitude(B_true) * np.array([ B_body[0], B_body[1], B_body[2] ])

    
    def propagatePosition(self,  controlForce, disturbanceForce, dt):
        self.t += dt
        
        self.controlForce = controlForce
        self.disturbanceForce = disturbanceForce
        
        # new_orbit = self.orbit.propagate(dt << u.second)
        new_orbit = fpropagate(self.orbit, controlForce, self.m, self.A, self.t, dt)
        
        self.longitude, self.latitude, self.altitude = ECI_to_ECEF((new_orbit.r << u.meter).value, self.t)
        
        self.orbit = new_orbit
        

class satelliteState_orbit():
    def __init__(self, satellite_orbit, satellite_mass, satellite_area):
        self.m = satellite_mass
        self.A = satellite_area
        
        self.a = satellite_orbit[0] << u.m
        self.e = satellite_orbit[1] << u.one
        self.i = satellite_orbit[2] << u.deg
        self.O = satellite_orbit[3] << u.deg
        self.w = satellite_orbit[4] << u.deg
        self.nu = satellite_orbit[5] << u.deg
        
        self.orbit = Orbit.from_classical(Earth, self.a, self.e, self.i, self.O, self.w, self.nu)
        
        self.controlForce = [0, 0, 0]
        self.disturbanceForce = [0, 0, 0]
        
        self.latitude = 0
        self.longitude = 0
        self.altitude = 0
        
        self.t = 0
  
    def propagatePosition(self,  controlForce, disturbanceForce, dt):
        self.t += dt
        
        self.controlForce = controlForce
        self.disturbanceForce = disturbanceForce
        
        # new_orbit = self.orbit.propagate(dt << u.second)
        new_orbit = fpropagate(self.orbit, controlForce, self.m, self.A, self.t, dt)
        
        self.longitude, self.latitude, self.altitude = ECI_to_ECEF((new_orbit.r << u.meter).value, self.t)
        
        self.orbit = new_orbit    



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

def ECEF_to_ECI(lla, t):
    
    transformer = pyproj.Transformer.from_crs({"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
                                              {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'})
    
    x, y, z = transformer.transform(lla[1], lla[0], lla[2], radians=False)

    gamma = (360/(23.9345*3600)) * t
    gamma *= np.pi/180
    gamma *= -1
    
    x_ = x*np.cos(gamma) - y*np.sin(gamma) # meters
    y_ = x*np.sin(gamma) + y*np.cos(gamma) # meters
    z_ = z # meters
    
    return [x_, y_, z_]

    

def appendDataFrame(df, state, t, scheduler, error, rwheel):
    tempdf = pd.DataFrame({'type':['attitude'],
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
                           'latitude':[state.latitude], 'longitude':[state.longitude], 'altitude':[state.altitude],
                           'aT':[scheduler._p['aT']/10**3], 'eT':[scheduler._p['eT']], 'iT':[scheduler._p['iT']*180/np.pi],
                           'q_dot':[scheduler.q_dot]
                           })
    return pd.concat([df, tempdf])


def appendDataFrame_orbit(df, state, scheduler, t):
    tempdf = pd.DataFrame({'type':['orbit'],
                           'time':[t], 'a':[(state.orbit.a << u.km).value],
                           'e':[state.orbit.ecc.value], 'i':[(state.orbit.inc << u.deg).value],
                           'aT':[scheduler._p['aT']/10**3], 'eT':[scheduler._p['eT']], 'iT':[scheduler._p['iT']*180/np.pi],
                           'x':[(state.orbit.r[0] << u.km).value], 'y':[(state.orbit.r[1] << u.km).value], 'z':[(state.orbit.r[2] << u.km).value],
                           'vx':[(state.orbit.v[0] << u.km/u.second).value], 'vy':[(state.orbit.v[1] << u.km/u.second).value], 'vz':[(state.orbit.v[2] << u.km/u.second).value],
                           'latitude':[state.latitude], 'longitude':[state.longitude], 'altitude':[state.altitude],
                           'q_dot':[scheduler.err]})
    return pd.concat([df, tempdf])

def appendDataFrame_spec(df, state, t, error):
    tempdf = pd.DataFrame({'type':['attitude'],
                           'time':[t], 
                           'Mx':[state.controlTorque[0]], 'My':[state.controlTorque[1]], 'Mz':[state.controlTorque[2]],
                           'Fx':[state.controlForce[0]], 'Fy':[state.controlForce[1]], 'Fz':[state.controlForce[2]],
                           'wx':[state.angularRate[0]], 'wy':[state.angularRate[1]], 'wz':[state.angularRate[2]],
                           'q1':[state.quaternion[0]], 'q2':[state.quaternion[1]], 'q3':[state.quaternion[2]], 'q4':[state.quaternion[3]],
                           'q1e':[error[0]], 'q2e':[error[1]], 'q3e':[error[2]], 'q4e':[error[3]],
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
    img = plt.imread("_resources/Earth.jpg")
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
        v.append( qm.magnitude([x[i], y[i], z[i]]) )
    
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


def plot_orbit_transfer(df):
    
    fig, ax = plt.subplots(2, 2)
    
    # generate plots for user
    ax[0, 0].set_title('semimajor axis vs. time')
    ax[0, 0].plot(df['time'], df['aT']-6371.8, 'k-', label='target semimajor axis (agl)')
    ax[0, 0].plot(df['time'], df['a']-6371.8, label='semimajor axis (agl)')
    ax[0, 0].set_xlabel('time (s)')
    ax[0, 0].set_ylabel('semimajor axis (agl, km)')
    ax[0, 0].legend()
    ax[0, 0].grid()

    ax[0, 1].set_title('eccentricity vs. time')
    ax[0, 1].plot(df['time'], df['eT'], 'k-', label='target eccentricity')
    ax[0, 1].plot(df['time'], df['e'], label='eccentricity')
    ax[0, 1].set_xlabel('time (s)')
    ax[0, 1].set_ylabel('eccentricity')
    ax[0, 1].legend()
    ax[0, 1].grid()

    ax[1, 0].set_title('inclination vs. time')
    ax[1, 0].plot(df['time'], df['iT'], 'k-', label='target inclination (deg)')
    ax[1, 0].plot(df['time'], df['i'], label='inclination (deg)')
    ax[1, 0].set_xlabel('time (s)')
    ax[1, 0].set_ylabel('inclination (deg)')
    ax[1, 0].legend()
    ax[1, 0].grid()
    
    ax[1, 1].set_title('error vs. time')
    ax[1, 1].plot(df['time'], df['q_dot'], 'k-', label='error')
    ax[1, 1].set_xlabel('time (s)')
    ax[1, 1].set_ylabel('error (unitless)')
    ax[1, 1].legend()
    ax[1, 1].grid()


    plt.show()

def plot_spec(df):
    t = np.array( df['time'] )
    
    q1e = np.array( df['q1e'] )
    q2e = np.array( df['q2e'] )
    q3e = np.array( df['q3e'] )
    q4e = np.array( df['q4e'] )
    

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
    plt.plot(t, df['Mx'], label='control torque (x)')
    plt.plot(t, df['My'], label='control torque (y)')
    plt.plot(t, df['Mz'], label='control torque (z)')
    
    plt.grid()
    plt.legend()
    
    plt.xlabel('time (s)')
    plt.ylabel('control torque (N-m)')

    ################

    plt.figure(3)
    plt.plot(t, df['wx'], label='wx')
    plt.plot(t, df['wy'], label='wy')
    plt.plot(t, df['wz'], label='wz')
    
    plt.grid()
    plt.legend()
    
    plt.xlabel('time (s)')
    plt.ylabel('angular rate (rad/s)')

    ################

    plt.show()