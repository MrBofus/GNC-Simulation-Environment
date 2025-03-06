import numpy as np
import numba as nb

@nb.jit(nopython=True)
def magnitude(vec):
    mag = 0
    for v in vec:
        mag += v**2
    
    return np.sqrt(mag)


@nb.jit(nopython=True)
def normalize(vec):
    
    m = magnitude(vec)
    
    v_ = []
    for v in vec:
        v_.append(v/m)
    
    return v_


@nb.jit(nopython=True)
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


@nb.jit(nopython=True)
def conjugate(q):
    return [ -q[0], -q[1], -q[2], q[3] ]


@nb.jit(nopython=True)
def quaternionIntegral(quaternion, angularRate, dt):
    
    wx, wy, wz = angularRate[0], angularRate[1], angularRate[2]
    qx, qy, qz, qw = quaternion[0], quaternion[1], quaternion[2], quaternion[3]
    
    return normalize( [ (0.5*dt) * ( wx*qw - wy*qz + wz*qy + qx ),
                        (0.5*dt) * ( wx*qz + wy*qw - wz*qx + qy ),
                        (0.5*dt) * (-wx*qy + wy*qx + wz*qw + qz ),
                        (0.5*dt) * (-wx*qx - wy*qy - wz*qz + qw ) ] )


@nb.jit(nopython=True)
def quaternionMultiply(q1, q2):
 
    w1, x1, y1, z1 = q1[3], q1[0], q1[1], q1[2]
    w2, x2, y2, z2 = q2[3], q2[0], q2[1], q2[2]
    
    return normalize( [ w1*x2 + w2*x1 + y1*z2 - y2*z1,
                        w1*y2 + w2*y1 - x1*z2 + x2*z1,
                        w1*z2 + w2*z1 + x1*y2 - x2*y1,
                        w1*w2 - x1*x2 - y1*y2 - z1*z2] )




@nb.jit(nopython=True)
def quaternionDifferenceToAngularVelocity(q1, q2, dt):
    return (2 / dt) * np.array([ q1[0]*q2[1] - q1[1]*q2[0] - q1[2]*q2[3] + q1[3]*q2[2],
                                 q1[0]*q2[2] + q1[1]*q2[3] - q1[2]*q2[0] - q1[3]*q2[1],
                                 q1[0]*q2[3] - q1[1]*q2[2] + q1[2]*q2[1] - q1[3]*q2[0] ])


@nb.jit(nopython=True)
def axis_to_quaternion(axis, rotation):
    return normalize( np.array([ axis[0]*np.sin(rotation/2),
                                 axis[1]*np.sin(rotation/2),
                                 axis[2]*np.sin(rotation/2),
                                 np.cos(rotation/2) ]) )

@nb.jit(nopython=True)
def quaternion_to_axis(quaternion):
    theta = 2*np.arccos(quaternion[3])
    return normalize( np.array([quaternion[0], 
                                quaternion[1], 
                                quaternion[2]]) / np.sin(theta/2) )


@nb.jit(nopython=True)
def dcm_to_quaternion(dcm):
    q1 = np.sqrt( 0.25*abs( 1 + dcm[0][0] - dcm[1][1] - dcm[2][2] ) )
    q2 = np.sqrt( 0.25*abs( 1 - dcm[0][0] + dcm[1][1] - dcm[2][2] ) )
    q3 = np.sqrt( 0.25*abs( 1 - dcm[0][0] - dcm[1][1] + dcm[2][2] ) )
    q4 = np.sqrt( 0.25*abs( 1 + dcm[0][0] + dcm[1][1] + dcm[2][2] ) )

    if q1 == max([ q1, q2, q3, q4 ]):
        q2 = (1/(4*q1)) * ( dcm[0][1] + dcm[1][0] )
        q3 = (1/(4*q1)) * ( dcm[2][0] + dcm[0][2] )
        q4 = (1/(4*q1)) * ( dcm[1][2] - dcm[2][1] )

        return normalize( np.array([q1, q2, q3, q4]) )

    if q2 == max([ q1, q2, q3, q4 ]):
        q1 = (1/(4*q2)) * ( dcm[0][1] + dcm[1][0] )
        q3 = (1/(4*q2)) * ( dcm[1][2] + dcm[2][1] )
        q4 = (1/(4*q2)) * ( dcm[2][0] - dcm[0][2] )

        return normalize( np.array([q1, q2, q3, q4]) )

    if q3 == max([ q1, q2, q3, q4 ]):
        q1 = (1/(4*q3)) * ( dcm[2][0] + dcm[0][2] )
        q2 = (1/(4*q3)) * ( dcm[1][2] + dcm[2][1] )
        q4 = (1/(4*q3)) * ( dcm[2][0] - dcm[0][2] )

        return normalize( np.array([q1, q2, q3, q4]) )

    if q4 == max([ q1, q2, q3, q4 ]):
        q1 = (1/(4*q4)) * ( dcm[1][2] - dcm[2][1] )
        q2 = (1/(4*q4)) * ( dcm[2][0] - dcm[0][2] )
        q3 = (1/(4*q4)) * ( dcm[0][1] - dcm[1][0] )

        return normalize( np.array([q1, q2, q3, q4]) )


def quaternion2euler(quaternion):
    qw, qx, qy, qz = quaternion[3], quaternion[0], quaternion[1], quaternion[2]

    phi = np.arctan2( 2*(qw*qx + qy*qz), 1 - 2*(qx**2 + qy**2) )
    theta = -np.pi/2 + 2*np.arctan2( np.sqrt(1 + 2*(qw*qy - qx*qz)), np.sqrt(1 - 2*(qw*qy - qx*qz)) )
    psi = np.arctan2( 2*(qw*qz + qx*qy), 1 - 2*(qy**2 + qz**2) )

    return phi, theta, psi


@nb.jit(nopython=True)
def orthogonalize(v, u):
    v = np.array(normalize(v))
    u = np.array(normalize(u))

    proj_v_onto_u = ( np.dot(v, u)/np.dot(u, u) ) * u

    return normalize( v - proj_v_onto_u )


@nb.jit(nopython=True)
def az_el_from_basis_vectors(basis1, basis2, basis3):
    basis1 = normalize(basis1)
    basis2 = normalize(basis2)
    basis3 = normalize(basis3)

    dcm = np.array([basis1, basis2, basis3])

    quaternion = dcm_to_quaternion(dcm)

    x, y, z = quaternion[0], quaternion[1], quaternion[2]

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2( y , x)

    return np.degrees(theta), np.degrees(phi)


def return_c(B, e):
    return np.cross(-e, normalize(B))