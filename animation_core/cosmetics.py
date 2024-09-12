import numpy as np


def interpolate_between_colors(color1, color2, i):
    return i*np.array(color1) + (1-i)*np.array(color2)


def plot_with_gradient(axes, c1, c2,
                       x, y, z,
                       gap, linewidth, gradient, skip_step):
    '''
    for i in range(1, len(x)):
        colorcode = [ 1 - abs(t[i-1]/max(t)), 
                     1, 
                     1 ]
        axes.plot(x[i-1:i+1], y[i-1:i+1], z[i-1:i+1], c=colorcode)
    '''
    
    for i in range(1, len(x), gap):
        '''
        iterator = abs(t[i]/max(t)) ** 0.5
        if iterator > 1: iterator = 0.999999
        '''

        iterator = i/len(x)

        '''
        colorcode = interpolate_between_colors([0.9607, 0.9843, 0.1764], 
                                               [0.5725, 0.0000, 0.5843], iterator)
        '''
        '''
        colorcode = interpolate_between_colors([0.9999, 0.9999, 0.1764], 
                                               [0.4000, 0.0000, 0.4000], iterator)
        '''

        colorcode = interpolate_between_colors(c1, c2, iterator)

        axes.plot(x[i-1:i+gap][::skip_step], y[i-1:i+gap][::skip_step], z[i-1:i+gap][::skip_step], c=colorcode, linewidth=linewidth, alpha=gradient)


