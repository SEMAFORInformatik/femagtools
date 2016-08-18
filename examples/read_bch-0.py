
import femagtools
import matplotlib.pylab as pl
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

bch = femagtools.read_bchfile('TEST_001.BCH')

beta = bch.ldq['beta']
i1 = bch.ldq['i1']
torque = bch.ldq['torque']

fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title(bch.project)
ax.set_xlabel(u'I1 [A]')
ax.set_ylabel(u'Beta [°]')
ax.set_zlabel(u'Torque [Nm]')
ax.azim = 210
X, Y = np.meshgrid(i1, beta)
ax.plot_surface(X, Y, torque,
                rstride=1, cstride=1, cmap=cm.jet)
pl.show()

