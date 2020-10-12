import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colorbar as cb
import matplotlib.colors as colors
import matplotlib as mpl
import numpy as np
import yt
import p64_productlog.phase_tools as phase_tools
reload(phase_tools)
from importlib import reload
import pdb

plt.clf()
plt.close('all')

if 1:
    prefix='ca03' #something to note the simulations you used.
    directory = '/scratch1/dcollins/Paper49_EBQU/ca03_turb_weak'
    frame = 339 
if 1:
    prefix='ca02' #something to note the simulations you used.
    directory = '/scratch1/dcollins/Paper49_EBQU/ca02_turb'
    frame=100
plot_directory = "%s/PlotsToTransfer"%os.environ['HOME']



#This conditional tests to see if region has been defined.
#if not, define it.  This process is slow so don't repeat if you don't have to.
ds = yt.load('%s/DD%04d/data%04d'%(directory,frame,frame))
region = ds.all_data()

#set up bins.
Nbins=64
density_min = region['density'].min()
density_max = region['density'].max()
density_bins=np.logspace(np.log10(density_min),np.log10(density_max),Nbins)
velocity_magnitude_min = region['velocity_magnitude'].min()
velocity_magnitude_max = region['velocity_magnitude'].max()
velocity_magnitude_bins=np.logspace(np.log10(velocity_magnitude_min),np.log10(velocity_magnitude_max),Nbins/4)
bins={'density':density_bins,'velocity_magnitude':velocity_magnitude_bins}

#
# Produce joint and marginalized distributions.
# Formally these are not PDFs since we don't explicityly normalize,
# BUT since the total volume is 1, they are in fact PDFs.
#
if 'joint' not in dir():
    bin_fields = ['density','velocity_magnitude']
    joint =     yt.create_profile(region, bin_fields=bin_fields,fields=['cell_volume'], weight_field=None, override_bins=bins)
    prof_x  = yt.create_profile(region, bin_fields=['density'],fields=['cell_volume'], weight_field=None, override_bins=bins)
    prof_y  = yt.create_profile(region, bin_fields=['velocity_magnitude'],fields=['cell_volume'], weight_field=None, override_bins=bins)

    #
    # Plot things the way YT things they should be, for reference
    #
    p_joint = yt.PhasePlot.from_profile(joint)
    p_joint.save("%s/%s"%(plot_directory,prefix))

    plot_mag = yt.ProfilePlot.from_profiles(prof_x)
    plot_mag.save("%s/%s"%(plot_directory,prefix))
    plot_theta = yt.ProfilePlot.from_profiles(prof_y)
    plot_theta.save("%s/%s"%(plot_directory,prefix))



#
# Check that P(V1,V2) = P(V1) P(V2)
#
if 1:

    fig3,ax3_list = plt.subplots(1,2)
    ax30=ax3_list[0]
    ax31=ax3_list[1]
    phase_tools.plot_phase_contours_1(fig3,ax30,joint,ax2=ax31)
    fig3.savefig("%s/%s_pf2_%04d.pdf"%(plot_directory,prefix,frame))

