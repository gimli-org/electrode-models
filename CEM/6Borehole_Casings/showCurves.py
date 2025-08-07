import numpy as np
import pygimli as pg
import matplotlib.pyplot as plt

H0 = pg.core.DataMap('hom0.collect')
H1 = pg.core.DataMap('hom1.collect')
I0 = pg.core.DataMap('inhom0.collect')
I1 = pg.core.DataMap('inhom1.collect')

data = pg.DataContainerERT( 'cpp.shm' )

rhoa0 = I0.data(data) / H0.data(data)
rhoa1 = I1.data(data) / H1.data(data)

x = np.array( [data.sensorPosition(i).x() for i in range(1,20)] )

fig, ax = plt.subplots()
ax.plot( x, rhoa0, x, rhoa1 )
ax.legend( ('Node','CEM'), loc='best' )
ax.grid()
ax.set_xlabel( 'pole-dipole separation [m]' )
ax.set_ylabel( r'$\rho_a$ in $\Omega$m' )
plt.show()