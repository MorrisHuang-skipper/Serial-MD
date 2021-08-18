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
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['font.size'] = 18
plt.rcParams['figure.figsize'] = [5.6*2, 4]
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 8
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1

colors = cm.get_cmap('Set1', 10)

fig = plt.figure(constrained_layout=True)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)

ax = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])

# ax = fig.add_subplot(1, 3, 1, projection='3d')

ax.xaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', right='on')
ax2.xaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', top='on')
ax2.xaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', top='on')
ax2.yaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', right='on')
ax2.yaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', right='on')

fname = 'ext'
rx = np.loadtxt('../../data/'+fname+'/rx.dat')
# ry = np.loadtxt('../../data/'+fname+'/ry.dat')
# rz = np.loadtxt('../data/'+fname+'/rz.dat')
# vx = np.loadtxt('../data/'+fname+'/vx.dat')
# vy = np.loadtxt('../data/'+fname+'/vy.dat')
# vz = np.loadtxt('../data/'+fname+'/vz.dat')
t, Ki, Ke, K, U, Te, Ti, within = np.loadtxt('../../data/'+fname+'/info.dat', unpack=True)

step = rx.shape[0]
NP = rx.shape[1]
e = 1.60217e-19
kt = e/38.9
nth = 250


ax.plot(t*1e15, (Ki+Ke)/e/NP, '-', label=r'$K_{avg}$', color=colors(2))
ax.plot(t*1e15, U/e/NP, '-', label=r'$U_{avg}$', color=colors(3))
ax.plot(t*1e15, Ki/e/NP+Ke/e/NP+U/e/NP, '-', label=r'$E_{avg}$', color=colors(4))
ax.set_xlabel('$time \ [fs]$')
ax.set_ylabel('$E \ [ev]$')
ax.set_xlim(0, 300)
ax.axvspan(0, 100, alpha=0.3, color=colors(1))
ax.set_title('energy')

ax2.plot(t*1e15, Te/e, '-', label=r'$T_e$', color=colors(6))
ax2.plot(t*1e15, Ti/e, '-', label=r'$T_i$', color=colors(0))
ax2.set_xlabel('$time \ [fs]$')
ax2.set_ylabel('$Temperature \ [eV]$')
#  ax3.set_ylabel('$Temperature \ [k_BT]$')
ax2.set_xlim(0, 300)
ax2.axvspan(0, 100, alpha=0.3, color=colors(1))
ax2.set_title('temperature')

ax.legend(loc=1)
ax2.legend(loc=1)

plt.tight_layout()
# plt.savefig('ext.png', dpi=2000)
plt.show()


# os.system('ffmpeg -i ../figures/'+fname+'/%04d.png -c:v ffv1 -qscale:v 0 ../figures/animate.mp4')

#  plt.tight_layout()
