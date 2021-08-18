import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pylab import cm
import math
from mpl_toolkits.mplot3d import Axes3D
import os
import sys
import matplotlib.gridspec as gridspec

mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['font.size'] = 16
plt.rcParams['figure.figsize'] = [5.6, 4]
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 6
plt.rcParams['legend.fontsize'] = 13
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1

colors = cm.get_cmap('Set1', 10)

fig = plt.figure(constrained_layout=True)

ax = fig.add_subplot(1, 1, 1)

# ax = fig.add_subplot(1, 3, 1, projection='3d')

ax.xaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', right='on')

fname = 'ext'
rx = np.loadtxt('../data/'+fname+'/rx.dat')
ry = np.loadtxt('../data/'+fname+'/ry.dat')
rz = np.loadtxt('../data/'+fname+'/rz.dat')
vx = np.loadtxt('../data/'+fname+'/vx.dat')
vy = np.loadtxt('../data/'+fname+'/vy.dat')
vz = np.loadtxt('../data/'+fname+'/vz.dat')
t, Ki, Ke, K, U, Te, Ti, within = np.loadtxt('../data/'+fname+'/info.dat', unpack=True)


step = rx.shape[0]
NP = rx.shape[1]
e = 1.60217e-19
kt = e/38.9
nth = 30

ecm = []
icm = []

# plot phase space
# for i in range(0, step, 1):
#     num = '{0:04}'.format(i)
    # ecm.append(ry[i, :-(NP//2)]*1e9/(NP//2)-20)
    # icm.append(ry[i, (NP//2):]*1e9/(NP//2)-20)

ax.set_xlim(7, 33)
ax.set_ylim(7, 33)
for i in range(50):
    ax.plot(rx[:, i]*1e9, ry[:, i]*1e9, '-')

plt.tight_layout()
plt.show()

