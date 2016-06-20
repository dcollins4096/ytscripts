

#
#Unit Tool.  To sort out prefactors of scaling.
#Pick density, temperature, tell me what G is.
# 
if 'oober' not in dir():
    GravitationalConstant = 1142 #p27
    GravitationalConstant = 607 #d122
else:
    GravitationalConstant = oober.pf.parameters['GravitationalConstant']

#GravitationalConstant = 1.62e3 #b02
#kelvin
T = 10
#in cm^{-3}
n0_cgs = 1000
density_marker = "%d"%n0_cgs

G_code = GravitationalConstant/(4*na.pi)
rho_code = 1.0
L_code = 1.0
t_code = 1.0
#amu/particle
mu = 2.3
amu = 1.660468e-24
rho_cgs = n0_cgs*mu*amu #grams/cc
G_cgs = 6.67e-8 #cm^3 g^-1 s^-2
GramsPerSolarMass = 1.98e33
#tff_cgs = na.sqrt(3*na.pi)
tff_cgs =  na.sqrt(3*na.pi/(32*G_cgs*rho_cgs) )
tff_cgs_myr = tff_cgs/(na.pi*1e7*1e6)
tff_code = na.sqrt( 3*na.pi/(32*G_code*rho_code) )
kb = 1.38065e-16 #gram cm^2/(second Kelvin)
cs_cgs = na.sqrt( T*kb/(amu*mu)) 
Tback = cs_cgs**2/kb*amu*mu
L_cgs = tff_cgs*cs_cgs/tff_code
pc_per_cm = 3.24077929e-19
cm_per_pc = 1./3.24077929e-19
L_pc = L_cgs * pc_per_cm
AvPerColumn = 1./1.9e21
AvBar = L_cgs*n0_cgs * AvPerColumn
print "L parsec %0.3e, n0 (cm^-3} %0.3e, T (K) %0.3e AvBar %0.2e"%(L_pc,n0_cgs,T,AvBar)


from yt.utilities import physical_constants as units
co_12_frequency = 115.271e9
mass_12 = (12.+16.)*units.mass_hydrogen_cgs
co_13_frequency = 110.201e9
mass_13 = (13.+16.)*units.mass_hydrogen_cgs
bandwidth = 25e6 #FCRAO, from Ridge et al 2006
n_bins = 1024    #FCRAO, from Ridge et al 2006

line_frequency = co_12_frequency
mass = mass_12

c_light = units.speed_of_light_cgs
kbolz = units.boltzmann_constant_cgs
mproton = units.mass_hydrogen_cgs
Distance = 150 * cm_per_pc
