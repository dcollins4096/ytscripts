import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import yt
#yt.enable_plugins()
if 1:
    not_in_yt3 = False #a universal flag to indicate things that need to be ported.
    execfile('go_lite')
#dummy
#Stuff I wrote that I like
#import davecallback
    import dsdiff_helpers
    import fPickle
    #execfile('xtra_fields.py')
    #execfile('xtra_energy_fields.py')
    execfile('davetools.py') #as needed, no more execfile.  
    import taxi
    if  not_in_yt3:
      import better_cbar
      #GET from img_region import * 
      #execfile('xtra_fields.py')
      execfile('scatter_fit.py')
      from uber import *
      #import stalker
      from time_marker import time_marker
      import tube
      import bilog
      import uberdiff

    from yt.analysis_modules.level_sets.api import * #for clumps
    import pyximport; pyximport.install()
#import particle_ops
    import random
#import clump_particles
#from yt.utilities.data_point_utilities import FindBindingEnergy
    from yt.utilities.physical_constants import \
                gravitational_constant_cgs as G
    import taxi

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
    MultiSpecies1=[ "Electron_Density", "HII_Density",\
        "HI_Density", "HeIII_Density", "HeII_Density", "HeI_Density"]
    "Cooling_Time",
    MultiSpecies2 = [ "Electron_Density", "H2II_Density", "H2I_Density", "HII_Density", "HI_Density",
        "HM_Density", "HeIII_Density", "HeII_Density", "HeI_Density"]

#ef('dave_callbacks.py')

