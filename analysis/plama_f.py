import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pylab import cm
import math
from mpl_toolkits.mplot3d import Axes3D
import os
import sys
import matplotlib.gridspec as gridspec
from scipy.fft import fft, fftfreq

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
spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)

ax = fig.add_subplot(spec[0, 0])

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
# ry = np.loadtxt('../data/'+fname+'/ry.dat')
# rz = np.loadtxt('../data/'+fname+'/rz.dat')
# vx = np.loadtxt('../data/'+fname+'/vx.dat')
# vy = np.loadtxt('../data/'+fname+'/vy.dat')
# vz = np.loadtxt('../data/'+fname+'/vz.dat')
t, Ki, Ke, K, U, Te, Ti, within = np.loadtxt('../data/'+fname+'/info.dat', unpack=True)


step = rx.shape[0]
n = rx.shape[1]
e = 1.60217e-19
kt = e/38.9
nth = 250
me = 9.10938e-31 
eps0 = 8.854187e-12



# rx = np.array(rx)
Lx = 40*1e-9
dt = 1e-15 / 10
tmp = 0
cnt = 0
# periodic condition
# while (abs(rx.sum()-tmp) > Lx/2):
#     print(abs(rx.sum()-tmp))
#     cnt += 1
#     print(cnt)
#     tmp = rx.sum()
# for i in range(0, n):
#     for j in range(1, step):
#         while (abs(rx[j, i]-rx[j-1, i]) >= Lx/2):
#             if ((rx[j, i] - rx[j-1, i])>0):
#                 rx[j, i] -= Lx 
#             else:
#                 rx[j, i] += Lx
# for i in range(0, n):
#     for j in range(1, step):
#         while (abs(ry[j, i]-ry[j-1, i]) >= Lx/2):
#             if ((ry[j, i] - ry[j-1, i])>0):
#                 ry[j, i] -= Lx 
#             else:
#                 ry[j, i] += Lx
# for i in range(0, n):
#     for j in range(1, step):
#         while (abs(rz[j, i]-rz[j-1, i]) >= Lx/2):
#             if ((rz[j, i] - rz[j-1, i])>0):
#                 rz[j, i] -= Lx 
#             else:
#                 rz[j, i] += Lx


ne = n / (40*1e-9)**3
omegap = (ne * e**2 / (me * eps0))**0.5
fp = omegap / (2 * np.pi)
print('plasma frequency', '{:e}'.format(fp))
print('period', 1/fp)

ax.plot(t, within, '-')
# ax.plot(t, 50*np.sin(2*np.pi*fp*t)+300)
# xf = fftfreq(step, dt)[:step//2]
# yf = fft(within)
# yf[0] = 0
# ax.plot(xf, 2/step*abs(yf[0:step//2]), '-.')
# ax.set_xlim(0, 3*fp)
# ax.vlines(x=fp, ymin=0, ymax=0.001, color='grey', alpha=0.5, ls='--')

# t = 0
# time = []
# xmean = []
# ymean = []
# zmean = []
# for i in range(0, step):
#     time.append(t)
#     xmean.append(sum(rx[i, :(n//2)])/n)
#     ymean.append(sum(ry[i, :(n//2)])/n)
#     zmean.append(sum(rz[i, :(n//2)])/n)
#     t += dt
# ax.plot(time, xmean, '.-')
# ax.plot(time, ymean, '.-')
# ax.plot(time, zmean, '.-')

# window = np.hamming(step)
# yf = fft(xmean)
# xf = fftfreq(step, dt)[:step//2]
# yf[0] = 0
# ax.plot(xf, 2/step*abs(yf[0:step//2]), '*-')

# yf = fft(ymean)
# xf = fftfreq(step, dt)[:step//2]
# yf[0] = 0
# ax.plot(xf, 2/step*abs(yf[0:step//2]), '*-')

# yf = fft(zmean)
# xf = fftfreq(step, dt)[:step//2]
# yf[0] = 0
# ax.plot(xf, 2/step*abs(yf[0:step//2]), '*-')

# ax.vlines(x=fp, ymin=0, ymax=4e-10, color='grey', alpha=0.5, ls='--')
# ax.set_xlim(0, 5*fp)


# xf = fftfreq(step, dt)[:step//2]
# t = np.arange(0, 50*1e-15, dt)
# for i in range(100):
#     ax.plot(time, rx[:, i], label=i)
#     # ax.plot(t, 1e-10*np.sin(2*np.pi*fp*t))
    # yf = fft(rx[:, i])
    # yf = fft(1e-10*np.sin(2*np.pi*fp*t))
    # ax.plot(xf, 2/step*abs(yf[0:step//2]), '-')
# ax.vlines(x=fp, ymin=0, ymax=2e-8, color='grey', alpha=0.5, ls='--')
# ax.set_xlim(-fp, 3*fp)


plt.tight_layout()
plt.show()
