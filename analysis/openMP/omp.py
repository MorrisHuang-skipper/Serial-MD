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

ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(4))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(5))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(2))

N = np.array([1, 2, 4, 8, 16, 32, 64, 128])
trd10 = np.array([72, 206, 713, 2360, 8731, 35453, 145688, 589580])
trd20 = np.array([114, 320, 956, 3140, 11530, 45623, 187895, 784535])
trd40 = np.array([260, 643, 1582, 4693, 15752, 60393, 251171, 1014813])

trd10 = trd10 / 10
trd20 = trd20 / 20
trd40 = trd40 / 40

# ax.plot(N, tlink, '*', color=colors(0))
# ax.plot(N, tmove, '*', color=colors(1))
# ax.plot(N, tout, '*', color=colors(2))
ax.plot(N, trd10/3600, '-*', color=colors(0), markersize=10, label='10 threads')
ax.plot(N, trd20/3600, '-s', color=colors(1), markersize=10, label='20 threads')
ax.plot(N, trd40/3600, '-^', color=colors(3), markersize=10, label='40 threads')
ax.set_xlabel(r'Number of particles ($\times 10^3$)', labelpad=10)
ax.set_ylabel(r'$CPU \ time \ [hr]$', labelpad=10)
ax.set_title('Computation time (threads number)')
ax.legend()


plt.tight_layout()
# plt.savefig('threads_time.eps')
# plt.savefig('threads_time.png', dpi=2000)
plt.show()
