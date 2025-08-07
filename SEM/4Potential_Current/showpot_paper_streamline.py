# -*- coding: utf-8 -*-

"""
Created on Wed May 23 14:12:36 2012
"""

import numpy as np
import pygimli as pg
from pygimli.mplviewer import drawField, drawStreams
from math import atan2
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, Normalize


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


#mesh = g.Mesh( 'BHmesh/mesh' )
#pot = g.RMatrix( 'BHpot/num.bmat' )
fig, ax = plt.subplots(4, 2, figsize=(6.9, 6))
fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
nPOT = [1, 2, 3]
maxZ = 10
for ii in nPOT:
    if ii < 3:
        mesh = pg.load('BHmesh2parts/mesh.bms')
        pot = pg.load('BHpot2parts/num'+str(ii)+'.bmat')
        diffPot = pot[0] - pot[2]
    else:
        mesh = pg.load('BHmesh3parts/mesh.bms')
        pot = pg.load('BHpot3parts/num'+str(ii)+'.bmat')
        diffPot = pot[0] - pot[3]

    print(mesh)
    print(mesh.dimension())
    mesh.createNeighbourInfos()
    grid2d= pg.createGrid(x=np.linspace(-250.,250., 50),
                          y=np.linspace(-200.,  0., 50))

    diffPot2d = pg.interpolate(mesh, diffPot,
                               x=pg.x(grid2d.positions()),
                               y=pg.x(grid2d.positions())*0.0,
                               z=pg.y(grid2d.positions()))

    diffPot2d /= abs(max(diffPot2d))
    diffPot2d *= 10.
    print(grid2d)
    pg.show(grid2d, diffPot2d, axes=ax[ii-1, 0], nLevs=21)
    #pg.show(mesh, diffPot*1000, axes=ax[ii-1,0])
    #drawField(ax[ii-1,0],mesh,diffPot*1000,nLevs=21,colors='black')

    maxZ = max(pg.abs(diffPot2d))
    ax[ii-1, 0].set_xlim([-250., 250.])
    ax[ii-1, 0].set_ylim([-200., 0.])

    name = 'rho.map'
    m, rho = np.loadtxt(name, comments='#', delimiter=' ', unpack=True)
    xBH = 150.
    zBH = 100.
    x_m = np.arange(-250., 250.5, 2)
    z_m = np.arange(0., -200.5, -2)
    D_abs = np.zeros(shape=(len(z_m), len(x_m)))
    D_phi = np.zeros(shape=(len(z_m), len(x_m)))
    D_x = np.zeros(shape=(len(z_m), len(x_m)))
    D_z = np.zeros(shape=(len(z_m), len(x_m)))
    for zz in range(len(z_m)):
        for xx in range(len(x_m)):
            pos = pg.RVector3(x_m[xx], 0., z_m[zz])
            c = mesh.findCell(pos, False)
            if c is not None:
                d = -((rho[c.marker()-1])**-1) * c.grad(pos, diffPot)
                D_x[zz, xx] = d[0]
                D_z[zz, xx] = d[2]
                D_abs[zz, xx] = np.sqrt(sum(d**2))
                D_phi[zz, xx] = atan2(d[2], d[0])

    mp.rcParams['font.size'] = 10
    extent = (x_m[0], x_m[len(x_m)-1], z_m[len(z_m)-1], z_m[0])
    X, Z = np.meshgrid(x_m, -z_m)
    inc, inc2 = 5, 10
    X1 = X[::inc, ::inc2]
    Z1 = Z[::inc, ::inc2]
    D_x1 = D_x/D_abs
    D_z1 = D_z/D_abs
    ax[ii-1, 0].tick_params('both', length=5, width=1.2, which='major', pad=7)

#    ax[ii-1,1].imshow( P.log(D_abs), vmin=-18,vmax=-10,
#                      interpolation="nearest", extent=extent )

    gridI2d = pg.createGrid(x=x_m, y=z_m)
    currentDens = np.asarray(np.log(D_abs+1e-8).flat)
    currentDens /= max(currentDens)
    currentDens = 1.0 - currentDens
    pg.show(gridI2d, currentDens, axes=ax[ii-1, 1], nLevs=256, omitLines=1)
    # streams start here .. für mehr streamlines das mesh hier feiner machen
    gridCoarse2d = pg.createGrid(x=np.linspace(-250., 250., 20),
                                 y=np.linspace(-200., 0., 20))

    print("coarseMesh for streams", gridCoarse2d)
    drawStreams(ax[ii-1, 1], grid2d, diffPot2d, coarseMesh=gridCoarse2d,
                color='black')
    if ii == 3:
        ax[ii-1, 0].set_xlabel('x [m]')
        ax[ii-1, 1].set_xlabel('x [m]')

#    ax[ii-1, 1].set_yticklabels([])
    plt.sca(ax[ii-1, 1])
    plt.yticks([-200, -150, -100, -50, -0], ['']*5)
    ax[ii-1, 1].tick_params('both', length=5, width=1.2, which='major', pad=7)

norm = Normalize(vmin=-maxZ, vmax=maxZ)
cbar = ColorbarBase(ax[3, 0], norm=norm, cmap=cm.jet, orientation='horizontal')
cbar.set_label('$u$ [mV]', fontsize=10)
cbar.ax.tick_params(labelsize=10)

norm = Normalize(vmin=-18, vmax=-10)
cbar = ColorbarBase(ax[3, 1], norm=norm, cmap=cm.jet, orientation='horizontal')
cbar.set_label(r'log(|$j$|) [A/m$^{2}$]', fontsize=10)
cbar.ax.tick_params(labelsize=10)

for aa in ax[3, :]:
    aa.set_aspect(0.07)
    po = np.array(aa.get_position().bounds)
    po[1] += 0.03
    aa.set_position(po)

plt.draw()

xsh, ysh = -315, 25
ax[0, 0].text(xsh, ysh, '(a)', fontsize=8)
ax[0, 1].text(xsh, ysh, '(b)', fontsize=8)
ax[1, 0].text(xsh, ysh, '(c)', fontsize=8)
ax[1, 1].text(xsh, ysh, '(d)', fontsize=8)
ax[2, 0].text(xsh, ysh, '(e)', fontsize=8)
ax[2, 1].text(xsh, ysh, '(f)', fontsize=8)

#fig.savefig('fig4.pdf', bbox_inches='tight')
#fig.savefig('fig4.eps', bbox_inches='tight')
fig.savefig('fig4.tif', bbox_inches='tight', dpi=200)

#pg.wait()
plt.show()
