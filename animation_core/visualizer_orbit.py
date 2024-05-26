#```````````````````````````````````````````````````````````````````````````````````````````````````````````````##
# written by: ME !! :)
#

from stl import mesh
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import math
import numpy as np
import copy
import animation_core.renderLighting as rl
from animation_core.cosmetics import plot_with_gradient
import sys
# sys.path.append("../GNC-Simulation-Environment")
import quaternion_math.quaternionMath as qm


# execute the mesh update on the plot to animate stl motion
# --> each time this runs it's like one "frame" of animation
def handleUpdate(df, ax, ax_a, ax_e, ax_i, index: int, scale: int):
    # clear current plot of any graphics
    ax.clear()
    ax_a.clear()
    ax_e.clear()
    ax_i.clear()
    
    # update title (example, not needed)
    ax.set_title('t = ' + str(round(df['time'].iloc[index]/(3600*24), 1)) + ' days')


    # ax.plot(df['x'].tail(index), df['y'].tail(index), df['z'].tail(index))
    plot_with_gradient(ax, [0.2, 0.2, 1.0], [0.1, 0.1, 0.5], 
                       df['x'].iloc[:index], df['y'].iloc[:index], df['z'].iloc[:index],
                       int(index/10 + 1), 1, 1)
    
    
    ammount = 20000
    if index < ammount:
        ax_a.set_title('semimajor axis (km)')
        ax_a.plot(df['time'].iloc[:index], df['a'].iloc[:index]-6378.1)
        ax_a.plot(df['time'].iloc[:index], df['aT'].iloc[:index]-6378.1, 'k-')
        ax_a.grid()
        
        ax_e.set_title('eccentricity')
        ax_e.plot(df['time'].iloc[:index], df['e'].iloc[:index])
        ax_e.plot(df['time'].iloc[:index], df['eT'].iloc[:index], 'k-')
        ax_e.grid()
        
        ax_i.set_title('inclination (deg)')
        ax_i.plot(df['time'].iloc[:index], df['i'].iloc[:index])
        ax_i.plot(df['time'].iloc[:index], df['iT'].iloc[:index], 'k-')
        ax_i.grid()
    
    else:
        ax_a.set_title('semimajor axis (km)')
        ax_a.plot(df['time'].iloc[(index-ammount):index], df['a'].iloc[(index-ammount):index]-6378.1)
        ax_a.plot(df['time'].iloc[(index-ammount):index], df['aT'].iloc[(index-ammount):index]-6378.1, 'k-')
        ax_a.grid()
        
        ax_e.set_title('eccentricity')
        ax_e.plot(df['time'].iloc[(index-ammount):index], df['e'].iloc[(index-ammount):index])
        ax_e.plot(df['time'].iloc[(index-ammount):index], df['eT'].iloc[(index-ammount):index], 'k-')
        ax_e.grid()
        
        ax_i.set_title('inclination (deg)')
        ax_i.plot(df['time'].iloc[(index-ammount):index], df['i'].iloc[(index-ammount):index])
        ax_i.plot(df['time'].iloc[(index-ammount):index], df['iT'].iloc[(index-ammount):index], 'k-')
        ax_i.grid()
    
    ax_a.set_ylim([min(df['a'])-6378.1-100, max(df['a'])-6378.1+100])
    ax_e.set_ylim([min(df['e'])-0.00001, max(df['e'])+0.00001])
    ax_i.set_ylim([min(df['i'])-0.5, max(df['i'])+0.5])
    
    ax_a.xaxis.set_ticklabels([])
    ax_e.xaxis.set_ticklabels([])
    
    # scale the graph
    ax.auto_scale_xyz([-scale, scale], 
                      [-scale, scale], 
                      [-scale, scale])
    
    
    
    ax.set_axis_off()
    ax.set_facecolor('black') 


def runVisualizer_orbit(df, pause_amount: float, **kwargs):
    if 'buff_amt' in kwargs.keys():
        buffer_amount = kwargs['buff_amt']
    else:
        buffer_amount = 1

    if 'scale' in kwargs.keys():
        scale = kwargs['scale']
    else:
        scale = 100

    # Create a new plot
    figure = plt.figure()
    # figure.set_facecolor('black')
    ax_3d = figure.add_subplot(3, 3, 1, projection='3d', position=[0, 0.15, 0.7, 0.7])
    ax_3d.view_init(elev=-10, azim=-120, roll=180)
    
    ax_altitude = figure.add_subplot(3, 3, 3)
    ax_eccentricity = figure.add_subplot(3, 3, 6)
    ax_inclination = figure.add_subplot(3, 3, 9)
    
    ax_altitude.set_title('altitude (km)')
    ax_eccentricity.set_title('eccentricity')
    ax_inclination.set_title('inclination (deg)')

    for t in range(len(df)):
        if t%buffer_amount == 0:
            handleUpdate(df, ax_3d, ax_altitude, ax_eccentricity, ax_inclination, t, scale)
            plt.pause(pause_amount)
    
    plt.show()


if __name__ == "__main__":
    import pandas as pd

    df = pd.read_csv( '_out/simulation_results.txt' )

    runVisualizer_orbit(df, 0.01, buff_amt=100)
