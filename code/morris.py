import f90nml
from math import floor, sqrt, pi

nano = 1.e-9
ang = 1.e-10
femto = 1.e-15
atto = 1.e-18
c0 = 2.99792458e8
eps0 = 8.854187e-12

#######################################################
# simulation dimension
dimension = 3
# radius of spherical target
radius = 10*nano
# number of particle of species one
np1 = 2000
# number of particle of species two
np2 = 1000
# simulation size
lx = 40*nano
ly = 40*nano
lz = 40*nano
# output location
loc = '../../data/ext/'
# timestep
dt = atto
# maximum simulation time
maxt = 100*femto
# particle width (Mic)
w0 = ang
# PIC level particle size
wpic = 5*ang
# cutoff radius (x3 wpic)
rcut = 3*wpic
# number of output times
dump = 100
# initial temperature
temperature = 0
# initial charge distribution (sphere, box)
init_distribution = 'sphere'
# initial velocity distribution (uniform, gaussian)
init_v_distribution = 'gaussian'
# periodic boundary for particles
periodic = False
#######################################################
# laser(external field) parameter
laser_switch = False
laser_freq = 1 / femto
laser_wavelength = c0 / laser_freq
laser_Intensity = 1.e16
laser_Emax = sqrt(2*laser_Intensity/c0/eps0)
laser_start = 50*femto
laser_end = 100*femto
#######################################################


#######################################################
nml = f90nml.read('InputFile/default.nml')

nml['input']['dim'] = dimension
nml['input']['radius'] = radius
nml['input']['np1'] = np1
nml['input']['np2'] = np2
nml['input']['lx'] = lx
nml['input']['ly'] = ly
nml['input']['lz'] = lz
nml['input']['loc'] = loc
nml['input']['dt'] = dt
nml['input']['maxt'] = maxt
nml['input']['w0'] = w0
nml['input']['wpic'] = wpic
nml['input']['rcut'] = rcut
nml['input']['numdump'] = dump
nml['input']['temperature'] = temperature
nml['input']['init_dis'] = init_distribution
nml['input']['init_v'] = init_v_distribution
nml['input']['periodic'] = periodic

# laser parameter
nml['laser']['e_ext'] = laser_switch
nml['laser']['freq'] = laser_freq
nml['laser']['omega'] = 2*pi*laser_freq
nml['laser']['lambda'] = laser_wavelength
nml['laser']['waven'] = 2*pi/laser_wavelength
nml['laser']['intensity'] = laser_Intensity
nml['laser']['emax'] = laser_Emax
nml['laser']['laser_start'] = laser_start
nml['laser']['laser_end'] = laser_end

# derived parameter
nml['derived']['np'] = np1+np2
nml['derived']['nstep'] = int(maxt/dt)
nml['derived']['dumpstep'] = int(int(maxt/dt)/dump)
nml['derived']['dumptime'] = int(int(maxt/dt)/dump)*dt
nml['derived']['ncellx'] = floor(lx/rcut)
nml['derived']['ncelly'] = floor(ly/rcut)
nml['derived']['ncellz'] = floor(lz/rcut)
nml['derived']['totcell'] = floor(lx/rcut)*floor(ly/rcut)*floor(lz/rcut)
nml['derived']['clx'] = lx/floor(lx/rcut)
nml['derived']['cly'] = ly/floor(ly/rcut)
nml['derived']['clz'] = lz/floor(lz/rcut)

nml.write('InputFile/input.nml', force=True)
#######################################################
