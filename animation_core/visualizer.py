#```````````````````````````````````````````````````````````````````````````````````````````````````````````````##
# written by: ME !! :)
#

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
import math
import numpy as np
import copy
import renderLighting as rl

# rotate stl mesh given a quaternion rotation
def rotateGivenQuaternion(np_mesh, quaternion):
    '''
    # check for singularity
    if quaternion[3] > 0.999:
        print('its me im here :))')
        quaternion[3] = 0.999
    '''
    # calculate angle and unit vector from quaternion
    # https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/index.htm
    angle = 2*math.acos(quaternion[3])
    x_ = quaternion[0] / math.sqrt(1 - quaternion[3]*quaternion[3])
    y_ = quaternion[1] / math.sqrt(1 - quaternion[3]*quaternion[3])
    z_ = quaternion[2] / math.sqrt(1 - quaternion[3]*quaternion[3])

    # update mesh given rotation
    np_mesh.rotate( np.array([ x_, y_, z_ ]), angle )

    return np_mesh



# execute the mesh update on the plot to animate stl motion
# --> each time this runs it's like one "frame" of animation
def handleUpdate(mesh_, translation, vvec, ax, quaternion: np.array, index: int, scale: int):
    # clear current plot of any graphics
    ax.clear()
    
    # update title (example, not needed)
    ax.set_title('t = ' + str(index/50))

    # execute quaternion rotation
    tempmesh = copy.deepcopy(mesh_.meshlist['deployed_mesh'])
    tempmesh = rotateGivenQuaternion(tempmesh, quaternion)


    tempmesh.translate(translation)
    deployed_vecs = tempmesh.vectors

    earth_vecs =    100*mesh_.meshlist['earth_mesh'].vectors


    # render the stl as a matplotlib object
    deployed_obj = mplot3d.art3d.Poly3DCollection(deployed_vecs, zorder=100)
    earth_obj =    mplot3d.art3d.Poly3DCollection(earth_vecs, zorder=0, alpha=1.0)



    # color the object grey and make the edges visible
    deployed_obj.set_edgecolor('k')
    deployed_obj.set_facecolor((0.7, 0.7, 0.7))

    rgbMesh = rl.renderLighting(mesh_.meshlist['earth_mesh'])
    earth_obj.set_facecolor(rgbMesh)

    ax.quiver( translation[0], translation[1], translation[2],
               vvec[0], vvec[1], vvec[2],
               length=25 )

    # add the render to the matplotlib figure
    ax.add_collection3d( deployed_obj )
    ax.add_collection3d( earth_obj )


    # scale the graph
    ax.auto_scale_xyz([translation[0]-scale, translation[0]+scale], 
                      [translation[1]-scale, translation[1]+scale], 
                      [translation[2]-scale, translation[2]+scale])
    
    ax.set_axis_off()
    ax.set_facecolor('black') 


class mesh_organizer_():
    def __init__(self, **kwargs):
        self.meshlist = {}
        for k in kwargs.keys():
           self.meshlist.update({k: kwargs[k]})



def runVisualizer(qlist: list, rlist:list, eventlist: list, pause_amount: float, **kwargs):
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
    axes = figure.add_subplot(projection='3d')
    axes.view_init(elev=-10, azim=-92, roll=0)

    # Load the STL files and add the vectors to the plot
    sat_mesh = mesh.Mesh.from_file('_resources/satellite_model_deployed.stl')
    earth_mesh = mesh.Mesh.from_file('_resources/earth_model.stl')

    meshhandler = mesh_organizer_(deployed_mesh=sat_mesh, 
                                  earth_mesh = earth_mesh)

    for t in range(len(qlist)):
        if t%buffer_amount == 0:
            handleUpdate(meshhandler, rlist[t], eventlist[t], axes, qlist[t], t, scale)
            plt.pause(pause_amount)


def temp(df):
    figure = plt.figure()
    axes = figure.add_subplot(projection='3d')

    for i in range(len(df)):
        axes.clear()

        axes.quiver( 0, 0, 0,
                     df['wx'].iloc[i], df['wy'].iloc[i], df['wz'].iloc[i],
                     length=25 )
        axes.quiver( 0, 0, 0,
                     df['Mx'].iloc[i], df['My'].iloc[i], df['Mz'].iloc[i],
                     length=25, color='r' )
        
        axes.auto_scale_xyz([-5, 5], 
                            [-5, 5], 
                            [-5, 5])
        
        plt.pause(0.01)


if __name__ == "__main__":
    import pandas as pd

    df = pd.read_csv( '_out/simulation_results.txt' )
    
    qlist = []
    rlist = []
    vlist = []
    for i in range(len(df)):
        qlist.append(np.array([ df['q1'].iloc[i], df['q2'].iloc[i], 
                                df['q3'].iloc[i], df['q4'].iloc[i] ]))
        
        rlist.append(np.array([ df['x'].iloc[i], df['y'].iloc[i], df['z'].iloc[i] ]))
        vlist.append(np.array([ df['vx'].iloc[i], df['vy'].iloc[i], df['vz'].iloc[i] ]))

    runVisualizer(qlist, rlist, vlist, 0.01, buff_amt=10)
