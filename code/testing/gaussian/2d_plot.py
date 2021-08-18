import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pylab import cm

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

colors = cm.get_cmap('Set1', 9)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.xaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', right='on')


e = 1.6e-19
x = np.loadtxt('out.dat', unpack=True)
ax.hist(x, color=colors(0), bins=500, histtype='step', density=True)
x = np.loadtxt('out2.dat', unpack=True)
ax.hist(x, color=colors(1), bins=500, histtype='step', density=True)
x = np.loadtxt('out3.dat', unpack=True)
ax.hist(x, color=colors(2), bins=500, histtype='step', density=True)

plt.tight_layout()
# plt.savefig('../figure/1a.pdf')
plt.show()
