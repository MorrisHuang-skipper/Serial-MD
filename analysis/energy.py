import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pylab import cm
import math

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

colors = cm.get_cmap('Set1', 5)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#  ax2 = fig.add_subplot(3, 1, 2)
#  ax3 = fig.add_subplot(3, 1, 3)

ax.xaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', right='on')

#  ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
#  ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
#  ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(.05))
#  ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(.01))

e = 1.60217e-19
epsilon = 8.854187e-12
Me = 9.10938356e-31
Mi = 1836 * Me
timestep = [r'$150\times 10^{-18}$', r'$100\times 10^{-18}$', r'$50 \times 10^{-18}$', r'$30 \times 10^{-18}$', '$0.1$', '$.2$']
for i in range(1, 5):
    fname = 'conv'+str(i)
    t, Ki, Ke, K, U, T, coll, accum, px, py, pz = np.loadtxt('../data/'+fname+'/info.dat', unpack=True)
    E = K + U
    #  ax.plot(t, Ki+Ke, label=r'$E$')
    #  ax.plot(t, U, label=r'$U$')
    #  ax.plot(t, E, label=r'$dt=$'+timestep[i], color=colors(i))
    #  ax.plot(t, px)
    #  ax.plot(t, py)
    #  ax.plot(t, pz)
    ax.plot(t*1e15, abs((E-E[0])/E[0])/coll, label=r'$dt=$'+timestep[i-1], color=colors(i))
    #  ax.plot(t, abs((E-E[0])/E[0]), label=r'$dt=$'+timestep[i], color=colors(i))
    #  ax2.plot(t, coll, '--', color=colors(i))
    #  ax3.plot(t, T, color=colors(i))
    #  ax2.plot(t, accum, '--', color=colors(i))

ax.set_yscale('log')
ax.set_xlabel(r'$Time \ [fs]$')
ax.set_ylabel(r'$\dfrac{\Delta E_{loss}}{collision} \ [\varepsilon/coll.]$')
ax.set_title('Mic3 ($R=1\AA, \ NP=1000, \, r_{cut}=3w_{pic}$)')

ax.legend()
plt.tight_layout()
plt.savefig('../figures/conv.eps')
plt.savefig('../figures/conv.png', dpi=1200)
plt.show()
