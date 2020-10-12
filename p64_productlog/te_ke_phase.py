

from go import *
import xtra_energy as xe
reload(xe)
import xtra_operators as xo
reload(xe)
import p64_productlog.phase_tools as phase_tools
reload(phase_tools)
car = taxi.load('ca07'); frame = 2100; sig=5
#car = taxi.load('ca05'); frame = 800; sig=2.3


def symlogspace(tw_min,tw_max,nbins,linthresh=1,linbins=2):
    thresh = np.log10(linthresh)
    minlog = np.log10( np.abs(tw_min))
    maxlog = np.log10( np.abs(tw_max))
    binsize = (maxlog-thresh)/nbins + (minlog-thresh)/nbins 
    out1 = -10**(np.arange( minlog,thresh,-binsize ))
    out2 = np.linspace(-linthresh,linthresh,linbins)
    out3 = 10**np.arange(np.log10(tw_max), thresh, -binsize)[::-1]
    return np.concatenate([out1,out2,out3])

def te_ke_ratio(field,data):
    return data['therm_energy']/data['kinetic_energy']

def v2(field,data):
    return 0.5*data['velocity_magnitude']**2
def therm_work(field,data):
    return xo.AdotDel(data,['velocity_%s'%s for s in 'xyz'], 'density')
if 'always' not in dir():
    always = False
car.derived_fields['acoustic'] = xe.add_energies
#nbins1 = nbins2 = 65
def plot_phase(ax,phase, field='cell_volume'):
    x_bins = phase.x_bins
    y_bins = phase.y_bins
    x_cen = 0.5*(x_bins[1:]+x_bins[:-1])
    y_cen = 0.5*(y_bins[1:]+y_bins[:-1])
    nbins1 = x_cen.size
    nbins2 = y_cen.size
    KIN = np.r_[(nbins2)*[x_cen]].transpose()
    TTT = np.r_[(nbins1)*[y_cen]]# .transpose()
    CCC = phase[field] #.transpose()  #transpose because pcolormesh plotting

    norm = mpl.colors.LogNorm(vmin=CCC[CCC>0].min(),vmax=CCC.max())
    cmap = mpl.cm.jet
    cmap.set_under('w')
    color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

    ploot=ax.pcolormesh(KIN,TTT,CCC,norm=norm)
    ax.set_xlabel(phase.x_field)
    ax.set_ylabel(phase.y_field)
    fig.colorbar(ploot,ax=ax)

if 'ds' not in dir() or always:
    car.frames=[frame]
    #directory = "/scratch1/dcollins/Paper49_EBQU/ca07_nofield_mach10"; frame=2100
    #ds = yt.load("%s/DD%04d/data%04d"%(directory, frame, frame))
    ds = car.load(frame)
    xe.add_energies(ds)
    ds.add_field( 'v2',function=v2,units='code_velocity**2')
    ds.add_field('therm_work',function=therm_work,units='code_density*code_velocity/code_length',
                validators=[yt.ValidateSpatial(1,'velocity-%s'%s) for s in 'xyz'],sampling_type='cell')
    ds.add_field('te_ke_ratio',function=te_ke_ratio,units='dimensionless',sampling_type='cell')
    region = ds.all_data()


if 1:

    #region['v2']
    therm = region['therm_energy']
    therm_min = therm.min()
    therm_max = therm.max()
    global_min = -1/np.e
    if therm_min < global_min:
        print("serious error in planning")
        raise


    nbins1 = 125
    nbins_lim=10
    #therm_bins = np.linspace(therm_min, therm_max, nbins1)
    if therm_max > np.abs(global_min):
        therm_bins1 = np.logspace( np.log10(np.abs(global_min)), np.log10(therm_max), nbins1)
        dx = (therm_bins1[1]-therm_bins1[0])
        dx = min( [dx, np.abs(2*therm_min/100)])
    else:
        therm_bins1=[]
        dx=(therm_max-therm_min)/nbins1
    if hasattr(dx,'v'):
        dx=dx.v
    #therm_bins2 = np.linspace(therm_min, np.abs(therm_min), nbins_lim).v
    therm_bins2 = np.arange( therm_min.v, np.abs(therm_min).v, dx)
    therm_bins = np.concatenate( [therm_bins2,therm_bins1])

    ke = region['kinetic_energy']
    ke_min = ke.min()
    ke_max = ke.max()
    nbins2 = 65
    #ke_bins = np.linspace(ke_min, ke_max, nbins2)
    ke_bins = np.logspace(np.log10(ke_min), np.log10(ke_max), nbins2)


if 'prof_therm_ke_dv' not in dir() or always:
    bins = {'therm_energy':therm_bins,'kinetic_energy':ke_bins}
    prof_therm_ke_dv  =  yt.create_profile(region, bin_fields=['therm_energy', 'kinetic_energy'],
                                        fields=['cell_volume'], weight_field=None, override_bins=bins)
    prof_therm_dv  =  yt.create_profile(region, bin_fields=['therm_energy'],
                                        fields=['cell_volume'], weight_field=None, override_bins=bins)


if 1:
    # 
    # TE KE phase diagram.
    # 
    fig,ax = plt.subplots(1,1)
    plot_phase(ax,prof_therm_ke_dv)
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xscale('symlog',linthreshx=0.3)

    outname = '../PigPen/%s_n%04d_phase_density_kinetic_energy_cell_volume'%(car.outname,frame)
    fig.savefig(outname)
    print(outname)
    plt.close(fig)

if 0:
    # 
    # Marginalized TE
    # 
    fig, axes = plt.subplots(1,2)
    ax0=axes[0]; ax1=axes[1]
    phase_tools.plot_phase_contours_1(fig,ax0,prof_therm_ke_dv, ax2=ax1)
    fig.savefig('../PigPen/%s_n%04d_contours.png'%(car.outname,frame))

    
