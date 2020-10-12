
from go import *
import xtra_energy as xe
reload(xe)
from scipy.special import lambertw  
import p64_productlog.ProductLog as PL
reload(PL)
#car = taxi.load('ca04'); frames=[200]; sig=0.1
#car = taxi.load('ca03_turb_weak'); frames = [300]; sig=2.3
car = taxi.load('ca05'); frames = [800]; sig=2.3
if 'always' not in dir():
    always = False
car.derived_fields['acoustic'] = xe.add_energies

def velocity_squared(field,data):
    return 0.5*data['velocity_magnitude']**2
def lrange(*args,**kwargs):
    return list(range(*args,**kwargs))
frames = [2,10,20]+lrange(50,350,50)
frames=[800]
always=True
for frame in frames:
    if 'prof_therm' not in dir() or always:
        car.frames=[frame]
        ds = car.load()
        ds.add_field( 'v2',function=velocity_squared,units='code_velocity**2',sampling_type='cell')
        region = ds.all_data()
        region['v2']
        min_dv = region['cell_volume'].min().v

        therm = region['therm_energy']
        therm_min = therm.min()
        therm_max = therm.max()
        bin_size = (therm_max-therm_min)/1000
        nbins = 1000
        therm_bins = np.linspace(therm_min, therm_max, nbins)
        bins = {'therm_energy':therm_bins}
        prof_therm  =  yt.create_profile(region, bin_fields=['therm_energy'],fields=['cell_volume'], weight_field=None, override_bins=bins)

    if 'prof_ke' not in dir() or always:
        ke = region['kinetic_energy']
        ke_min = ke.min()
        ke_max = ke.max()
        bin_size = (ke_max-ke_min)/1000
        nbins = 1000
        #ke_bins = np.linspace(ke_min, ke_max, nbins)
        bins = {'kinetic_energy': therm_bins} #ke_bins}
        prof_ke  =  yt.create_profile(region, bin_fields=['kinetic_energy'],fields=['cell_volume'], weight_field=None, override_bins=bins)

    if 'prof_v2' not in dir() or always:
        v2 = region['v2']
        v2_min = v2.min()
        v2_max = v2.max()
        nbins = 1000
        v2_bins = np.linspace(v2_min, v2_max, nbins)
        bins = {'v2': v2_bins} #v2_bins}

        prof_v2  =  yt.create_profile(region, bin_fields=['v2'],fields=['cell_volume'], weight_field=None, override_bins=bins)

    if 'prof_vel' not in dir() or always:
        vel = region['velocity_magnitude']
        vel_min = vel.min()
        vel_max = vel.max()
        nbins = 1000
        vel_bins = np.linspace(-vel_max, vel_max, nbins)
        bins={'velocity_magnitude':vel_bins}
        prof_vel  =  yt.create_profile(region, bin_fields=['velocity_magnitude'],fields=['cell_volume'], weight_field=None, override_bins=bins)
        vel = region['velocity_magnitude']
        dv = region['cell_volume']
        volume = dv.sum()
        vrms = np.sqrt((vel**2 * dv).sum()/volume)
    if 'prof_vx' not in dir() or always:
        vel = region['velocity_magnitude']
        vel_min = vel.min()
        vel_max = vel.max()
        nbins = 1000
        vel_bins = np.linspace(-vel_max, vel_max, nbins)
        bins={'velocity_x':vel_bins}
        bins['velocity_y']=vel_bins
        bins['velocity_z']=vel_bins
        prof_vx  =  yt.create_profile(region, bin_fields=['velocity_x'],fields=['cell_volume'], weight_field=None, override_bins=bins)
        prof_vy  =  yt.create_profile(region, bin_fields=['velocity_y'],fields=['cell_volume'], weight_field=None, override_bins=bins)
        prof_vz  =  yt.create_profile(region, bin_fields=['velocity_z'],fields=['cell_volume'], weight_field=None, override_bins=bins)


#if 'prof_ge' not in dir() or always:
#    prof_ge   =  yt.create_profile(region, bin_fields=['grav_energy'],fields=['cell_volume'], weight_field=None)#, override_bins=bins)

    if 'prof_rho' not in dir() or always:
        prof_rho  =  yt.create_profile(region, bin_fields=['density'],fields=['cell_volume'], weight_field=None)#, override_bins=bins)
        rho = region['density']
        dv = region['cell_volume']

    if 1:
        volume = dv.sum()
        mu_rho = (rho*dv).sum()/volume
        sigma_rho = np.sqrt(( (rho-mu_rho)**2*dv).sum()/volume)
        mu_lnrho = (np.log(rho/mu_rho)*dv).sum()/volume
        lnrho  = np.log(rho/mu_rho)
        sigma_lnrho = np.sqrt( ( (lnrho - mu_lnrho)**2*dv).sum()/volume)

        

    therm_bins = prof_therm.x_bins
    therm_cen = 0.5*(therm_bins[1:]+therm_bins[:-1])

    plt.close('all')
    fig,ax=plt.subplots(2,2)
    ax0=ax[0][0];ax1=ax[0][1]
    ax2=ax[1][0];ax3=ax[1][1]
    ax0.plot(therm_cen,prof_therm['cell_volume'],c='k',label='data',marker='*')
#r'$P(\rho \epsilon)$'

    scale=1
    x=therm_bins.v
    MachNumber = vrms.v
    b=1.0
    sigma_lnrho_predict = np.sqrt( np.log( 1+b**2*MachNumber**2))
    mu_lnrho_predict = -0.5*sigma_lnrho**2    

    predict = PL.g1( x,mu=mu_lnrho, sigma=sigma_lnrho,c=1,rho0=1,rho1=1).real
#scale = prof_therm['cell_volume'].max()/predict.max()
#scale = 1./predict.sum()
    scale = prof_therm['cell_volume'][6]/predict[6]
    ax0.plot( x, predict*scale,c='g',marker='*',label='almost data')
    ax0.set_ylim(min_dv/2, 1)

    predict = PL.g1( x,mu=mu_lnrho_predict, sigma=sigma_lnrho_predict,c=1,rho0=1,rho1=1).real
    predict[ predict<0] = 0
#scale = prof_therm['cell_volume'].max()/predict.max()
#scale = 1./predict.sum()
    ax0.plot( x, predict*scale,c='r',marker='*',label='Mach number only')
    ax0.set_ylim(min_dv/2, 1)

    axbonk(ax0,xlabel=r'$\rho \epsilon$',ylabel=r'$P(\rho \epsilon)$', yscale='log')
    ax0.set_xscale('symlog',linthreshx = 1)
    ax0.legend(loc=0)

    ke_bins = prof_ke.x_bins
    ke_cen = 0.5*(ke_bins[1:]+ke_bins[:-1])

    v2_bins = prof_v2.x_bins
    v2_cen = 0.5*(v2_bins[1:]+v2_bins[:-1])


    scale = prof_ke['cell_volume'][300]/prof_therm['cell_volume'][300]
    ax1.plot(ke_cen, prof_ke['cell_volume'], c='r', label='KE')
    ax1.plot(therm_cen, prof_therm['cell_volume']*scale, c='y', label='AE')
    ax1.plot(v2_cen, prof_v2['cell_volume'], c='g', label='AE')
    axbonk(ax1,ylabel=r'$P(E)$',xlabel=r'E (KE,AE)', yscale='log', xscale='linear')

    fig2,ax5=plt.subplots(1,1)
    scale = prof_ke['cell_volume'][300]/prof_therm['cell_volume'][300]
    ax5.plot(ke_cen, prof_ke['cell_volume'], c='r', label='KE')
    ax5.plot(therm_cen, prof_therm['cell_volume']*scale, c='y', label='AE')
    ax5.plot(v2_cen, prof_v2['cell_volume'], c='g', label='AE')
    axbonk(ax5,ylabel=r'$P(E)$',xlabel=r'E (KE,AE)', yscale='log', xscale='linear')
    fig2.savefig('../PigPen/%s_n%04d_AE_KE_log.png'%(car.outname,frame))


    for nprof,prof in enumerate([prof_vx,prof_vy,prof_vz]):
        vx_bins = prof.x_bins
        vx_cen = 0.5*(vx_bins[1:]+vx_bins[:-1])
        ax2.plot(vx_cen, prof['cell_volume'],c=[0.5]*4)#, c='g', label='AE')
    vmag_bins = prof_vel.x_bins
    vmag_cen = 0.5*(vmag_bins[1:]+vmag_bins[:-1])
    ax2.plot(vmag_cen, prof_vel['cell_volume'],c='k')#, c='g', label='AE')
        
    gss = np.exp(-(vx_cen)**2/(2*sig**2))
    gss *= prof['cell_volume'].max()/gss.max()
    ax2.plot(vx_cen, gss,c='b')
    maxwellian = vx_cen**2*gss
    ok = vx_cen > 0
    maxwellian *= prof_vel['cell_volume'].max()/maxwellian.max()
    ax2.plot(vx_cen[ok], maxwellian[ok],c='r')


    rho_bins = prof_rho.x_bins
    rho_cen = 0.5*(rho_bins[1:]+rho_bins[:-1])
    ax3.plot(rho_cen,prof_rho['cell_volume'])
    axbonk(ax3,xlabel=r'$\rho$',xscale='log',yscale='log')

#axbonk(ax3,ylabel=r'$P(E)$',xlabel=r'E (KE,AE)', yscale='linear',xscale='linear')
#ax3.plot(ke_cen, prof_ke['cell_volume'], c='r', label='KE')
#ax3.plot(therm_cen, prof_therm['cell_volume'], c='y', label='AE')
#axbonk(ax3,ylabel=r'$P(E)$',xlabel=r'E (KE,AE)', yscale='log',xscale='log')
##ax3.set_xscale('symlog',linthreshx=0.1)

    fig.savefig('../PigPen/%s_n%04d_acoustic_PDFs.png'%(car.outname,frame))
