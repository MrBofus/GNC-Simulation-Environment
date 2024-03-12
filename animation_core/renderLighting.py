from matplotlib.colors import LightSource
import numpy as np

# https://stackoverflow.com/questions/56864378/how-to-light-and-shade-a-poly3dcollection
# https://stackoverflow.com/questions/54006078/selection-of-face-of-a-stl-by-face-normal-value-threshold



def renderLighting(np_mesh):
    # https://github.com/WoLpH/numpy-stl/issues/163

    # Create light source
    ls = LightSource(azdeg=90, altdeg=0)

    # Darkest shadowed surface, in rgba
    dk = np.array([0.1, 0.1, 0.2, 1.0])
    # Brightest lit surface, in rgba
    lt = np.array([0.3, 0.3, 0.9, 1.0])
    # Interpolate between the two, based on face normal
    shade = lambda s: (lt-dk) * s + dk

    # Set face colors 
    sns = ls.shade_normals(np_mesh.get_unit_normals(), fraction=1.0)
    rgba = np.array([shade(s) for s in sns])

    return rgba