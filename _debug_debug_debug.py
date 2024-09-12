import numpy as np
import quaternion_math.quaternionMath as qm
import matplotlib.pyplot as plt

figure = plt.figure()
axes = figure.add_subplot(projection='3d')

"""
a x b = c

b = (c x a) / magnitude(a)**2 + t*a


        c x B = e
       -(B x c) = e
        B x c = -e
==>
        c = (-e x B) / magnitude(B)**2 + t*B
        
        
"""



B = np.array([1., -2., 3.])
e = np.array([4., 5., 6.])

c = np.cross(-e, B) / np.dot(B, B)**2 + 0*B

print('\n', e, '\n', np.cross(c, B) + (0.866)*e, '\n')

'''
B = np.array([1, 2, 3])

e = np.array([0.2, -0.3, 0.4])

c = np.cross(-e, B) / qm.magnitude(B)**2 + 0*B
'''

'''
axes.quiver(0, 0, 0, 
            B[0], B[1], B[2], length=1, color='r')
axes.quiver(0, 0, 0, 
            e[0], e[1], e[2], length=10, color='g')
axes.quiver(0, 0, 0, 
            c[0], c[1], c[2], length=10, color='b')

axes.auto_scale_xyz([-5, 5], 
                    [-5, 5], 
                    [-5, 5])

plt.show()
'''

'''
figure = plt.figure()
axes = figure.add_subplot(projection='3d')

axes.quiver(0, 0, 0, 
            e[0], e[1], e[2], length=10, color='g')
axes.quiver(0, 0, 0, 
            np.cross(c, B)[0], np.cross(c, B)[1], np.cross(c, B)[2], length=10, color='b')

axes.auto_scale_xyz([-5, 5], 
                    [-5, 5], 
                    [-5, 5])

plt.show()
'''