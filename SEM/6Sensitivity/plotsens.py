# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 14:11:47 2013

@author: Ronczka.M
"""

import pygimli as pg
import pybert as pb
import numpy as np
from scipy.ndimage import zoom
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase


def logdrop(value, eps=1e-6):
    logdata = np.zeros(len(value))

    for ii in range(len(value)):
        if value[ii] >= eps:
            logdata[ii] = np.log(abs(value[ii])/eps)
        elif value[ii] <= -eps:
            logdata[ii] = -np.log(abs(value[ii])/eps)
        else:
            logdata[ii] = 0

    return logdata

# load mesh and potential from calculation
hA = [1, 100, 100, 1, 100]  # hight of 1. borehole
hB = [1, 100, 100, 100, 1]  # hight of 2. borehole
hM = [1, 100, 1, 100, 50]
hN = [1, 100, 1, 1, 1]

#cdict = {'red': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)),
#         'green': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
#         'blue': ((0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0))}
cols = [(0, 0, 0), (1, 1, 1), (0, 0, 0)]
cmap = mp.colors.LinearSegmentedColormap.from_list('my_colormap', N=255,
                                                   colors=cols)

plt.rcParams['font.size'] = 10
fig, ax = plt.subplots(3, 2, figsize=(6.9, 6))

x_m = np.arange(-250., 250.5, 2)
z_m = np.arange(0., -200.5, -2)
XM, ZM = np.meshgrid(x_m, z_m)
x_p = np.array(XM.flat)
z_p = np.array(ZM.flat)
y_p = pg.RVector(len(x_p))
extent = [x_m[0], x_m[len(x_m)-1], z_m[len(z_m)-1], z_m[0]]

for a in range(5):
    data = pb.DataContainerERT('data_DpDp'+str(a)+'.ohm')
    mesh = pg.Mesh('meshDpDp/mesh'+str(a)+'.bms')
    mesh.createNeighbourInfos()

    f = pb.DCMultiElectrodeModelling(mesh, data, True)
    f.region(1).setBackground(True)
    f.createRefinedForwardMesh(True, True)

    mesh_inv = f.regionManager().paraDomain()

    k = 1. / data.get('r')
    data.set('k', k)
    vol = mesh_inv.cellSizes()

    model = pg.RVector(f.regionManager().parameterCount(), 1.0)  # homogeneous

    f.createJacobian(model)
    J = f.jacobian()
    print J.cols(), J.rows(), len(vol)
    sens2 = J[0] / vol

    for ii in range(len(J)):
        zwerg = logdrop(J[ii])
        out_sens = pg.asvector(zwerg)

    zwerg2 = logdrop(sens2, 1e-8)
    out_sens2 = pg.asvector(zwerg2)

    sens_z2 = pg.interpolate(mesh_inv, out_sens2, x_p, y_p, z_p) / 10
    sens_p2 = zoom(np.reshape(sens_z2, XM.shape), 4)

# plotting ####################################################################

    A = -hA[a]
    B = -hB[a]
    M = -hM[a]
    N = -hN[a]
    cl = 'k'
    aa = ax.flat[a]
#    aa = ax[int(a/2.), a % 2]
    aa.imshow(sens_p2, interpolation='nearest', extent=extent,
              cmap=cmap, vmin=-1, vmax=1)
#    co = aa.contour(zoom(np.flipud(sens_p2), 5), extent=extent,
#                    levels=[-0.4, -0.2, 0.0, 0.2, 0.4], colors='black')
#    plt.clabel(co, inline=1, fontsize=8, fmt='%g')
    aa.tick_params('both', length=5, width=1.2, which='major', pad=7)

    aa.plot([-150, -150], [0, A], cl, [-100, -100], [0, B], cl,
            [100, 100], [0, M], cl, [150, 150], [0, N], cl, linewidth=1.2)
    aa.set_xlabel('x [m]')
    aa.set_xlim(-250, 250)
    aa.text(-150, 5, 'A', va='bottom', ha='center')
    aa.text(-100, 5, 'B', va='bottom', ha='center')
    aa.text(100, 5, 'M', va='bottom', ha='center')
    aa.text(150, 5, 'N', va='bottom', ha='center')
    aa.text(-200, -30, '_', fontsize=16, va='center', ha='center')
    aa.text(-125, -40, '+', fontsize=16, ha='center')
    aa.text(-5, -30, '_', fontsize=16, ha='center')
    aa.text(120, -40, '+', fontsize=16, ha='center')
    aa.text(200, -30, '_', fontsize=16, ha='center')
    aa.set_yticks([-200, -150, -100, -50, 0])
#    plt.sca(aa)
#    plt.yticks([-200, -150, -100, -50, 0],['','','','',''])
    if a % 2 == 0:
        aa.set_ylabel('z [m]')
    else:
        aa.set_yticklabels(['']*5)

norm = Normalize(vmin=-1, vmax=1)
cbar = ColorbarBase(ax[2, 1], norm=norm, cmap=cmap, orientation='horizontal',
                    ticks=[-1, 0, 1])
cbar.set_ticklabels(['-1', '0', '+1'])
cbar.set_label('Sensitivity', fontsize=10)
cbar.ax.tick_params(labelsize=10)
ax[2, 1].set_aspect(0.07)

ms = 10
m = 'vk'
# mk surface electrode points
ax[0, 0].plot([-150, -100, 100, 150], [0, 0, 0, 0], m, ms=ms, zorder=5)
ax[1, 0].plot([100, 150], [0, 0], m, ms=ms, zorder=5)
ax[1, 1].plot([-150, 150], [0, 0], m, ms=ms, zorder=5)
ax[2, 0].plot([-100, 150], [0, 0], m, ms=ms, zorder=5)

xsh = -330
ysh = 25
ax[0, 0].text(xsh, ysh, '(a)')
ax[0, 1].text(xsh, ysh, '(b)')
ax[1, 0].text(xsh, ysh, '(c)')
ax[1, 1].text(xsh, ysh, '(d)')
ax[2, 0].text(xsh, ysh, '(e)')
