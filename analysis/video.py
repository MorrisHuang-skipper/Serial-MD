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
plt.rcParams['figure.figsize'] = [5.4*4, 5.2]
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 6
plt.rcParams['legend.fontsize'] = 13
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1

colors = cm.get_cmap('Set1', 10)

fig = plt.figure(constrained_layout=True)
spec = gridspec.GridSpec(ncols=4, nrows=1, figure=fig)

ax = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[0, 2])
ax4 = fig.add_subplot(spec[0, 3])

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
rx = np.loadtxt('../../data/'+fname+'/rx.dat')
ry = np.loadtxt('../../data/'+fname+'/ry.dat')
rz = np.loadtxt('../../data/'+fname+'/rz.dat')
vx = np.loadtxt('../../data/'+fname+'/vx.dat')
vy = np.loadtxt('../../data/'+fname+'/vy.dat')
vz = np.loadtxt('../../data/'+fname+'/vz.dat')
t, Ki, Ke, K, U, Te, Ti, within = np.loadtxt('../../data/'+fname+'/info.dat', unpack=True)

x1 = np.arange(-10, 10, 0.00001)
yy1 = (40-x1**2)**0.5
yy2 = -(40-x1**2)**0.5

step = rx.shape[0]
NP = rx.shape[1]
e = 1.60217e-19
kt = e/38.9
nth = 50

# plot phase space
for i in range(0, step, 1):
    num = '{0:04}'.format(i)
    # ax.set_zlim(0, 10)
    # 2d
    ax.plot(rx[i, :-(NP//2)]*1e9, ry[i, :-(NP//2)]*1e9, '.', color=colors(0), label='$H^+$', markersize=3)
    ax.plot(rx[i, (NP//2):]*1e9, ry[i, (NP//2):]*1e9, '.', color=colors(1), label='$e^{-}$', markersize=3)
    ax.plot(rx[:i, nth]*1e9, ry[:i, nth]*1e9, '-', color=colors(2), label='tracking')
    # for j in range(50):
    #     ax.plot(rx[:i, nth+j]*1e9, ry[:i, nth+j]*1e9, '.')

    ax.set_xlabel('$x \ [nm]$')
    ax.set_ylabel('$y \ [nm]$')
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 40)

    # diag
    ax2.plot(t[:i]*1e15, (Ki[:i]+Ke[:i])/e/NP, '-', label=r'$K_{avg}$', color=colors(2))
    ax2.plot(t[:i]*1e15, U[:i]/e/NP, '-', label=r'$U_{avg}$', color=colors(3))
    ax2.plot(t[:i]*1e15, Ki[:i]/e/NP+Ke[:i]/e/NP+U[:i]/e/NP, '-', label=r'$E_{avg}$', color=colors(4))
    ax2.set_xlabel('$time \ [fs]$')
    # ax2.set_ylabel('$E \ [ev]$')
    ax2.set_xlim(0, 100)

    ax3.plot(t[:i]*1e15, Te[:i]/e, '-', label=r'$T_e$', color=colors(6))
    ax3.plot(t[:i]*1e15, Ti[:i]/e, '-', label=r'$T_i$', color=colors(0))
    ax3.set_xlabel('$time \ [fs]$')
    ax3.set_ylabel('$Temperature \ [eV]$')
    #  ax3.set_ylabel('$Temperature \ [k_BT]$')
    ax3.set_xlim(0, 100)

    ax4.hist(vx[i, 1:NP//2], histtype='step')
    ax4.hist(vx[i, NP//2+1:NP], histtype='step')
    # ax4.set_xlabel('$time \ [fs]$')
    # ax4.set_ylabel('$Temperature \ [eV]$')
    # ax4.set_xlim(0, 100)

    ax.legend(loc=1)
    ax2.legend(loc=1)
    ax3.legend(loc=1)
    ax4.legend(loc=1)

    plt.tight_layout()
    plt.show(block=False)
    plt.pause(.01)
    # plt.savefig('../figures/'+fname+'/'+str(num)+'.png', dpi=300)
    ax.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    


# os.system('ffmpeg -i ../figures/'+fname+'/%04d.png -c:v ffv1 -qscale:v 0 ../figures/animate.mp4')

#  plt.tight_layout()
