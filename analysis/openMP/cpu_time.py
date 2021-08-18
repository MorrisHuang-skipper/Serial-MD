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
ax2 = ax.twinx()

ax.xaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', right='on')
ax2.yaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', right='on')
ax2.yaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', right='on')

ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(5000))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1000))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(50))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10))

N = [100, 200, 800, 2000, 4000, 8000, 20000]
md = [0.64, 1.92, 24.6, 156, 575, 2372, 14197]
mic = [0.30, 0.49, 3.8, 17, 51, 179, 1147]

# ax.plot(N, tlink, '*', color=colors(0))
# ax.plot(N, tmove, '*', color=colors(1))
# ax.plot(N, tout, '*', color=colors(2))
ax.plot(N, np.array(md)/60, '-*', color=colors(2), markersize=10, label='$MD$')
ax.plot(N, np.array(mic)/60, '-*', color=colors(1), markersize=10, label='$Mic$')
ax.set_xlabel('Number of particles', labelpad=10)
ax.set_ylabel(r'$CPU \ time \ [min]$', labelpad=10)
ax.set_title('Computation time: Mic3 vs. MD3')
ax.legend()

ax2.plot(N, np.array(md)/3600*0.7, '-*', color=colors(2), markersize=10, label='$MD$')
ax2.plot(N, np.array(mic)/3600*0.7, '-*', color=colors(1), markersize=10, label='$Mic$')
ax2.set_ylabel('NTD')

plt.tight_layout()
plt.savefig('cpu_time.eps')
plt.savefig('cpu_time.png', dpi=2000)
plt.show()
