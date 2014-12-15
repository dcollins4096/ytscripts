import yt
not_in_yt3 = False #a universal flag to indicate things that need to be ported.
execfile('go_lite')
#dummy
#Stuff I wrote that I like
#import davecallback
if  not_in_yt3:
  import better_cbar
  #GET from img_region import * 
  execfile('davetools.py') #as needed, no more execfile.  
  execfile('xtra_fields.py')
  execfile('scatter_fit.py')
  from uber import *
  #import stalker
  from time_marker import time_marker
  import tube
  import bilog
  import uberdiff


isothermal_hydro=['Density']+['%s-velocity'%s for s in 'xyz'] 
isothermal_mhd = isothermal_hydro + ['B%s'%s for s in 'xyz']
if 0:
    magnetic_fields = ['MagneticField_F_%s'%d for d in '123'] 
    electric_fields = ['ElectricField_%s'%d for d in '123']
else:
    magnetic_fields = ['BxF', 'ByF', 'BzF']
    electric_fields = ['Ex','Ey','Ez']
staggered = magnetic_fields + electric_fields
adiabatic_hydro = isothermal_hydro + ['TotalEnergy']
adiabatic_mhd = isothermal_mhd + ['TotalEnergy']
all_fields = adiabatic_mhd+staggered
alphabet='abcdefghijklmnopqrstuvwxyz'
