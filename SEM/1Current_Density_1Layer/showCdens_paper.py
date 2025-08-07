# -*- coding: utf-8 -*-
"""
Created on Wed Oct 08 15:58:27 2014

@author: Ronczka.M
"""

import pygimli as pg
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mp
from pygimli.viewer import show
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pygimli.mplviewer import drawStreamLines2

###############################################################################
def pointDataToBoundaryData(mesh, u):
    ret = np.zeros(mesh.boundaryCount())
    for b in mesh.boundaries():
        for i in range(b.nodeCount()):
            ret[b.id()] += u[b.node(i).id()]
        ret[b.id()] /= b.nodeCount()
    return ret
###############################################################################
def showPot2D3layer(mesh, uPot, axes, rho, title):
    grid = pg.createGrid(x=np.linspace(-5, 5, 50), y=np.linspace(0, -1000, 200))
    #grid.createNeighbourInfos()
    u = pg.interpolate(mesh, uPot[0], x=pg.x(grid.positions()),
                                      y=pg.z(grid.positions()),
                                      z=pg.y(grid.positions()))

    #r = pg.interpolate(mesh, rho, x=pg.x(grid.positions()),
    #                              y=pg.z(grid.positions()),
    #                              z=pg.y(grid.positions()))
#    print min(r),max(r)
#    r[pg.find(r<1)]=10.
#    show(grid, u, axes=axes, showLater=1)
    drawStreamLines2(axes, grid, -u)
    E = grid.boundaryDataToCellGradient(pointDataToBoundaryData(grid, u))
    #R = grid.boundaryDataToCellGradient(pointDataToBoundaryData(grid, r))
    Ix = -pg.x(E) / 10. #pg.x(R)
    Iy = -pg.y(E) / 10. #pg.y(R)
    for c in grid.cells():
        if c.center().y() < -300. and c.center().y() > -600.:
            Ix[c.id()] = Ix[c.id()] * 3
            Iy[c.id()] = Iy[c.id()] * 3

    Iabs = pg.sqrt(Ix*Ix + Iy*Iy)

    probline=np.linspace(grid.ymin(), grid.ymax(), 100)
    uProbe = pg.interpolate(grid, Iabs, x=probline*0.0 + 2.5, y=probline)

    return probline, uProbe
###############################################################################
def showPot2D1layer(mesh, uPot, axes, title):
    grid = pg.createGrid(x=np.linspace(-5, 5, 50), y=np.linspace(0, -60, 90))
    #grid.createNeighbourInfos()
    u = pg.interpolate(mesh, uPot[0], x=pg.x(grid.positions()),
                                      y=pg.z(grid.positions()),
                                      z=pg.y(grid.positions()))

#    show(grid, u, axes=axes, showLater=1)
    #drawStreamLines2(axes, grid, -u)
    E = grid.boundaryDataToCellGradient(pointDataToBoundaryData(grid, u))
    Ix = -pg.x(E) / 100.0
    Iy = -pg.y(E) / 100.0

    Iabs = pg.sqrt(Ix*Ix + Iy*Iy)
    ax, cbar = show(grid, Iabs, cMin=0.001, cMax=0.05, axes=axes, showLater=1, colorBar=True, orientation='vertical')
    cbar.ax.tick_params(labelsize=12)

    axes.set_title("$\\vec{j}$ and $|\\vec{j}|$ " + title, fontsize=20)

    divider = make_axes_locatable(axes)
    ax2 = divider.append_axes("bottom", size=1., pad=0.5)

    probline=np.linspace(grid.ymin(), grid.ymax(), 100)
    uProbe = pg.interpolate(grid, Iabs, x=probline*0.0 + 1.0, y=probline)

    ax2.plot(probline, uProbe)

    ax2.tick_params(labelsize=8)
    ax2.set_xlim([grid.ymax(), grid.ymin()])
    ax2.set_xlabel('$z$ -coordinate in m', fontsize=16)
    ax2.set_ylabel('$|\\vec{j}|$ in $A/m^2$', fontsize=16)

    return probline, uProbe
###############################################################################

# 3 layer #####################################################################
name = 'rho.map'
mCEM, rhoCEM = np.loadtxt( 'CEM_3layer/'+name, comments='#', delimiter = ' ', unpack=True )

#fig = plt.figure()
#ax1 = fig.add_subplot(1,3,1)
#ax2 = fig.add_subplot(1,3,2)
#ax3 = fig.add_subplot(1,3,3)

meshCEM = pg.Mesh( 'CEM_3layer/BHmesh/mesh.bms' )
potCEM = pg.RMatrix( 'CEM_3layer/BHpot/num.bmat' )
rho = [rhoCEM[c.marker()-1] for c in meshCEM.cells()]
print len(rho),min(rho),max(rho)
#d3, uCEM3 = showPot2D3layer(meshCEM, potCEM, ax1, np.array(rho), title='CEM' )

meshCC  = pg.Mesh( 'CC_3layer/BHmesh/mesh.bms' )
potCC  = pg.RMatrix( 'CC_3layer/BHpot/num.bmat' )
#d3, uCC3 = showPot2D3layer(meshCC, potCC, ax2, np.array(rho), title='CCM')

meshLS  = pg.Mesh( 'SEM_3layer/BHmesh/mesh.bms' )
potLS  = pg.RMatrix( 'SEM_3layer/BHpot/num.bmat' )
#d3, uSEM3 = showPot2D3layer(meshLS, potLS, ax3, np.array(rho), title='Line')
#plt.close()
###############################################################################
# 1 layer #####################################################################
#fig = plt.figure()
#ax1 = fig.add_subplot(1,4,1)
#ax2 = fig.add_subplot(1,4,2)
#ax3 = fig.add_subplot(1,4,3)
#ax4 = fig.add_subplot(1,4,4)

meshCEM1 = pg.Mesh( 'CEM_1layer/BHmesh/mesh' )
potCEM1 = pg.RMatrix( 'CEM_1layer/BHpot/num.bmat' )
#d1, uCEM1 = showPot2D1layer(meshCEM1, potCEM1, ax1, title='CEM')

meshCC1  = pg.Mesh( 'CC_1layer/BHmesh/mesh' )
potCC1  = pg.RMatrix( 'CC_1layer/BHpot/num.bmat' )
#d1, uCC1 = showPot2D1layer(meshCC1, potCC1, ax3, title='CC')

meshLS1  = pg.Mesh( 'SEM_1layer/BHmesh/mesh' )
potLS1  = pg.RMatrix( 'SEM_1layer/BHpot/num.bmat' )
#d1, uSEM1 = showPot2D1layer(meshLS1, potLS1, ax4, title='Line')
#plt.close()
zana1,jana1 = np.loadtxt( 'cdens_anal_d0_2.txt',unpack=True )

###############################################################################

plt.rcParams['font.size'] = 10
fig, ax = plt.subplots( 2,1,figsize=(3.33,5.5) )
lw=1.
cCEM, cCCM, cSEM, cEll = 'lightgray', 'darkgray', 'black', 'black'
ax[0].plot(-d1, uCEM1,cCEM,label='CEM',linewidth=lw )
ax[0].plot(-d1, uCC1,cCCM,label='CCM',linewidth=lw, linestyle='dashed' )
ax[0].plot(-d1, uSEM1,cSEM,label='SEM',linewidth=lw, linestyle='dotted' )
ax[0].plot( zana1,jana1,cEll,label='ellipse',linewidth=lw )
plt.sca(ax[0])
plt.yticks(np.arange(0,5e-3,0.5e-3),('0','','1e-3','','2e-3','','3e-3','','4e-3','') )
ax[0].set_ylabel(r'j [A/m$^ {2}$]',labelpad=5,fontsize=10)
#ax[0].set_xlabel('depth [m]',labelpad=7,fontsize=10)
ax[0].legend(loc=0)
ax[0].grid(True)
ax[0].tick_params('both', length=5, width=1.2, which='major',pad=7)

ax[1].plot(-d3, uCEM3,'c',label='CEM',linewidth=lw )
ax[1].plot(-d3, uCC3,'m',label='CCM',linewidth=lw )
ax[1].plot(-d3, uSEM3,'g',label='SEM',linewidth=lw )
test = [0,0.5e-4,1e-4,1.5e-4,2e-4,2.5e-4,3e-4,3.5e-4,4e-4,4.5e-4]
plt.sca(ax[1])
plt.yticks(test,('0','','1e-4','','2e-4','','3e-4','','4e-4','') )
plt.ylabel(r'j [A/m$^ {2}$]',labelpad=5,fontsize=10)
plt.xlabel('Depth [m]',labelpad=5,fontsize=10)
ax[1].set_xlim(0,1150)
ax[1].legend(loc=0)
ax[1].grid(True)
ax[1].tick_params('both', length=5, width=1.2, which='major',pad=7)

ax[0].text(-2.5, 4.8e-3, '(a)', fontsize=10)
ax[1].text(-100, 4.8e-4, '(b)', fontsize=10)
fig.subplots_adjust( top=0.92,left=0.25,right=0.9,bottom=0.09,wspace=0.2,hspace=0.29 )
plt.savefig('cdens_alles.pdf')#,bbox_inches='tight')
plt.savefig('cdens_alles.eps')
plt.show()