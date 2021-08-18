import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pylab import cm
import math

mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['font.size'] = 16
plt.rcParams['figure.figsize'] = [5.6, 4]
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 6
plt.rcParams['legend.fontsize'] = 15
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

ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(.2))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(.2))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(.05))

def Vmic(r, w0, wpic):
    return 1/r * (math.erf(r/np.sqrt(2)/w0) - math.erf(r/np.sqrt(2)/wpic))
def V(r):
    return 1/r
def F(r):
    return 1/r**2
def Fmicpic(r, w0, wpic):
    return -(np.sqrt(2)*np.exp(-r**2/(2*w0**2))/(np.sqrt(np.pi)*w0*r) - math.erf(r/np.sqrt(2)/w0)/r**2)\
           +(np.sqrt(2)*np.exp(-r**2/(2*wpic**2))/(np.sqrt(np.pi)*wpic*r) - math.erf(r/np.sqrt(2)/wpic)/r**2)
def Fmic(r, w0, wpic):
    return -(np.sqrt(2)*np.exp(-r**2/(2*w0**2))/(np.sqrt(np.pi)*w0*r) - math.erf(r/np.sqrt(2)/w0)/r**2)
           #  +(np.sqrt(2)*np.exp(-r**2/(2*wpic**2))/(np.sqrt(np.pi)*wpic*r) - math.erf(r/np.sqrt(2)/wpic)/r**2)

x, y1, y2 = [], [], []
y3, y4, y5, y6 = [], [], [], []
w0 = 0.5
wpic = 1.15
for i in np.arange(0.001, 5, .0001):
    x.append(i)
    y1.append(Vmic(i, w0, wpic))
    y2.append(V(i))
    y3.append(F(i))
    y4.append(Fmic(i, w0, wpic))
    y5.append(Fmicpic(i, w0, wpic))
    y6.append(Fmic(i, w0, wpic)-Fmicpic(i, w0, wpic))

ax.vlines(3.45, 0, 1, ls='--', color='grey', lw=2.5)
ax.axvspan(3.45, 5, alpha=0.3, color='grey')

ax.plot(x, y3, label='Coulomb', color='black')
ax.plot(x, y4, label='Mic', color=colors(0))
ax.plot(x, y5, label='MicPIC', color=colors(1))
ax.plot(x, y6, label='PIC', color=colors(2))

ax.text(1.3, .68, 'Coulomb', color='black', fontsize=18)
ax.text(.1, .88, 'MicPIC', color=colors(0), fontsize=18)
ax.text(.45, .58, 'Mic', color=colors(1), fontsize=18)
ax.text(.5, .18, 'PIC', color=colors(2), fontsize=18)
ax.text(2.9, .90, 'MD', color='black', fontsize=18)
ax.text(3.6, .90, 'PIC', color='black', fontsize=18)
ax.text(3.6, .50, r'$r_{cut}$', color='black', fontsize=18)

ax.set_ylim(0, 1)
ax.set_xlim(0, 5)
ax.grid()

ax.set_ylabel(r'$Force \quad [e^2/4\pi\epsilon_0]$')
ax.set_xlabel(r'$Distance \quad [PIC \ cells]$')
#  ax.legend()

plt.tight_layout()
plt.savefig('MicPIC.png', dpi=3000)
plt.show()
