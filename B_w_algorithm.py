import numpy as np


def transpose(vec):
    vec_candidate = []
    for v in vec:
        vec_candidate.append([v])
    
    return np.array(vec_candidate)


def magnitude(vec):
    magnitude = 0
    for v in vec:
        magnitude += v**2
    
    return np.sqrt(magnitude)

def normalize(vec):
    m = magnitude(vec)
    
    vec_candidate = []
    for v in vec:
        vec_candidate.append(v/m)
    
    return vec_candidate


def hamiltonProduct(v1, v2):
    return np.array([ v1[0]*v2[1] + v1[1]*v2[0] + v1[2]*v2[3] - v1[3]*v2[2],
                      v1[0]*v2[2] - v1[1]*v2[3] + v1[2]*v2[0] + v1[3]*v2[1],
                      v1[0]*v2[3] + v1[1]*v2[2] - v1[2]*v2[1] + v1[3]*v2[0],
                      v1[0]*v2[0] - v1[1]*v2[1] - v1[2]*v2[2] - v1[2]*v2[2] ])



def conjugate(vec):
    vec_candidate = []
    for v in vec:
        vec_candidate.append(-1 * v)
    
    vec_candidate[-1] = -1*vec_candidate[-1]
    
    return vec_candidate


def dcm2quat(dcm):
    if 1 + dcm[0][0] - dcm[1][1] - dcm[2][2] > 0:
        q1 = np.sqrt( 0.25 * (1 + dcm[0][0] - dcm[1][1] - dcm[2][2]) )
    else:
        q1 = - np.sqrt( 0.25 * abs(1 + dcm[0][0] - dcm[1][1] - dcm[2][2]) )
    
    if 1 - dcm[0][0] + dcm[1][1] - dcm[2][2] > 0:
        q2 = np.sqrt( 0.25 * (1 - dcm[0][0] + dcm[1][1] - dcm[2][2]) )
    else:
        q2 = - np.sqrt( 0.25 * abs(1 - dcm[0][0] + dcm[1][1] - dcm[2][2]) )
    
    if 1 - dcm[0][0] - dcm[1][1] + dcm[2][2] > 0:
        q3 = np.sqrt( 0.25 * (1 - dcm[0][0] - dcm[1][1] + dcm[2][2]) )
    else:
        q3 = - np.sqrt( 0.25 * abs(1 - dcm[0][0] - dcm[1][1] + dcm[2][2]) )
    
    if 1 + dcm[0][0] + dcm[1][1] + dcm[2][2] > 0:
        q4 = np.sqrt( 0.25 * (1 + dcm[0][0] + dcm[1][1] + dcm[2][2]) )
    else:
        q4 = - np.sqrt( 0.25 * (1 + dcm[0][0] + dcm[1][1] + dcm[2][2]) )
    
    term_1_2 = dcm[0][1]
    term_2_1 = dcm[1][0]
    term_1_3 = dcm[0][2]
    term_3_1 = dcm[2][0]
    term_2_3 = dcm[1][2]
    term_3_2 = dcm[2][1]
    
    if q1 == max([abs(q1), abs(q2), abs(q3), abs(q4)]):
        q = np.array([ 4 * q1**2,
                       term_1_2 + term_2_1,
                       term_1_3 + term_3_1,
                       term_2_3 - term_3_2 ])
        
    elif q2 == max([abs(q1), abs(q2), abs(q3), abs(q4)]):
        q = np.array([ term_1_2 + term_2_1,
                       4 * q2**2,
                       term_2_3 + term_3_2,
                       term_3_1 - term_1_3 ])
        
    elif q3 == max([abs(q1), abs(q2), abs(q3), abs(q4)]):
        q = np.array([ term_3_1 + term_1_3,
                       term_2_3 + term_3_2,
                       4 * q3**2,
                       term_1_2 - term_2_1])
        
    elif q4 == max([abs(q1), abs(q2), abs(q3), abs(q4)]):
        q = np.array([ term_2_3 - term_3_2,
                       term_3_1 - term_1_3,
                       term_1_2 - term_2_1,
                       4 * q4**2])
    
    return normalize(q)


def return_quaternion_from_measurements(Bm, Bc):
    R = Bm * transpose(Bc)
    return R, dcm2quat(R)

Bm = normalize([1, 0, 0])
Bc = normalize([0, 0, 1])

Bm_ = np.array([Bm[0], Bm[1], Bm[2], 0])
Bc_ = np.array([Bc[0], Bc[1], Bc[2], 0])

R, quaternion = return_quaternion_from_measurements(Bm, Bc)

print(R)
print(quaternion)

print(hamiltonProduct(hamiltonProduct(quaternion, Bc_), conjugate(quaternion)))