from mpl_toolkits import mplot3d
import numpy as np
from matplotlib import pyplot as plt
from tribarmes import TriangleBarMesh

ax = None

def set_aspect_equal_3d():
    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)
    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])
    print(xlim, ylim, zlim, plot_radius)
    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])
    
def loadplottriangles(fname):
    global ax
    tbm = TriangleBarMesh(fname)
    vs = tbm.GetBarMeshTriangles()
    vv = np.array(vs).ravel()
    X, Y, Z = vv[0::3], vv[1::3], vv[2::3]
    tris = np.array(range(len(vs)*3))
    tris.resize((len(vs), 3))
    plt.figure(figsize=(9,9))
    ax = plt.gca(projection='3d', proj_type='ortho')
    #ax.set_aspect('equal')

    ax.plot_trisurf(X,Y,Z, triangles=tris, linewidth=None)
    return tbm

def plottoolshape(ballrad, shaftheight, px, py, pz):
    tr = np.linspace(0, ballrad+1, ballrad*10+1)
    theta = np.linspace(0, 2*np.pi, 40)
    tr, theta = np.meshgrid(tr, theta)
    def toolheight(tr):
        return ballrad - np.sqrt(np.maximum(0, ballrad**2 - tr**2)) + \
                         np.maximum(0, (tr-ballrad)*(shaftheight))
    X = np.minimum(tr, ballrad) * np.sin(theta)
    Y = np.minimum(tr, ballrad) * np.cos(theta)
    Z = toolheight(tr)
    #ax = plt.gca(projection='3d', proj_type = 'ortho')
    from matplotlib.tri import Triangulation
    tri = Triangulation(np.ravel(tr), np.ravel(theta))
    tris = tri.triangles
    ax.plot_trisurf(X.ravel()+px, Y.ravel()+py, Z.ravel()+pz, triangles=tris, linewidth=None, color="red")
    return ax
