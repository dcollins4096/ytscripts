from go import *
reload(taxi)
import multi_imshow
reload(multi_imshow)

if 'car' not in dir():
    #car = taxi.load('b02')
    car = taxi.load('ca03_turb_weak')
car.outname = "%s/%s"%('plots_to_sort',car.name)
def costhetaz(field,data):
    return data['magnetic_field_z']/data['magnetic_field_strength']
def thetaz(field,data):
    output = np.zeros_like( data['costhetaz'])
    output[:] = np.arccos(data['costhetaz'][:].v)
    return output


def add_cos_theta(obj):
    obj.add_field('costhetaz',costhetaz,units='dimensionless')
    obj.add_field('thetaz',thetaz,units='dimensionless')

car.frames = [339]
car.derived_fields['add_cos_theta']=add_cos_theta
car.load()
##car.load(100)
##car.ds.all_data()['density']
#g=car.ds.index.grids[0]
plt.clf()

if 'prof_theta' not in dir():
    #car.fields=['costhetaz']
    #car.plot()
    Nbins=64
    ybins=np.linspace(0,np.pi*2,2*Nbins)
    xbins=np.logspace(-3,2,Nbins)
    bins={'thetaz':ybins,'magnetic_field_strength':xbins,'magnetic_field_z':xbins}
    #car.profile(['thetaz','cell_volume'],scales=['linear','linear'])
    phase=car.phase(['magnetic_field_strength','thetaz','cell_volume'],phase_args={'override_bins':bins})
    car.last_phase_plot.set_log('thetaz',False)
    car.last_phase_plot.save(car.outname)
    car.profile(['thetaz','cell_volume'],scales=['linear','linear'],override_bins=bins)
    prof_theta=car.last_prof
    car.profile(['magnetic_field_strength','cell_volume'],scales=['linear','linear'],override_bins=bins)
    prof_b2=car.last_prof
    car.profile(['magnetic_field_z','cell_volume'],scales=['linear','linear'],override_bins=bins)
    prof_bz=car.last_prof

if 1:
    fig,ax=plt.subplots(1,3)
    ax0=ax[0];ax1=ax[1]; ax2=ax[2]
    ph = phase[0]['cell_volume']
    ax0.plot(prof_theta['cell_volume'],c='k')
    ax0.plot(np.sum(ph,axis=0),c='r')
    ax0.set_title('theta')

    ax1.plot(prof_b2['cell_volume'],'k')
    ax1.plot(np.sum(ph,axis=1),'r')
    ax1.set_title('B2')

    b_bins = prof_bz.x_bins
    b_cen = 0.5*(b_bins[1:]+b_bins[:-1])
    dbin = b_cen[1:]-b_cen[:-1]
    pbz=prof_bz['cell_volume']
    dpbz=(pbz[1:]-pbz[:-1])/dbin
    ax2.plot(b_cen, b_cen*pbz)
    ax2.plot(b_cen[1:], -b_cen[1:]*dpbz,c='r')
    ax2.plot(b_cen, prof_b2['cell_volume'],c='k')
    ax2.set_xscale('log')
    plt.savefig('plots_to_sort/theta_dists.png')
        
if 1:
    p_theta=prof_theta['cell_volume']
    p_b2 = prof_b2['cell_volume']
    #p_b2.shape = (p_b2.size,1)
    p_theta.shape = (p_theta.size,1)
    ph2 = p_theta*p_b2
    ph=np.transpose(phase[0]['cell_volume'])
    multi_imshow.plotter2([ph,ph2],'plots_to_sort/p63_joint.png',labels=['2d','Kludge'],npx=2,norm='positive')

if 0:
    """works dont touch"""
    #car.fields=['costhetaz']
    #car.plot()
    ybins=np.arange(0,np.pi*2,0.01)
    xbins=np.arange(1e-3,100)
    bins={'costhetaz':ybins,'magnetic_field_strength':xbins}
    #car.profile(['thetaz','cell_volume'],scales=['linear','linear'])
    bins={'thetaz':ybins,'magnetic_field_strength':xbins}
    phase=car.phase(['magnetic_field_strength','thetaz','cell_volume'])
    #,phase_args={'override_bins':bins})
    car.profile(['thetaz','cell_volume'],scales=['linear','linear'])
    p1=car.last_prof
    car.profile(['magnetic_field_strength','cell_volume'],scales=['linear','linear'])
    p2=car.last_prof

if 0:
    #car.fields=['costhetaz']
    #car.plot()
    ybins=np.arange(0,np.pi*2,0.01)
    xbins=np.logspace(-3,2,50)
    bins={'costhetaz':ybins,'magnetic_field_strength':xbins}
    bins={'thetaz':ybins,'magnetic_field_strength':xbins}
    ext={'thetaz':[1,2],'magntic_field_strength':[1e-3,100]}
    phase=car.phase(['magnetic_field_strength','thetaz','cell_volume'] \
              ,phase_args={'override_bins':bins,'extrema':ext})
    phase[0].set_log('thetaz',False)
    phase[0].save(car.outname)
    
    #car.profile(['costhetaz','cell_volume'],scales=['linear','linear'])

if 0:
    plt.clf()
    prof = car.last_prof
    bin_edges = prof.x_bins
    bins =0.5*(bin_edges[1:]+bin_edges[:-1])
    plt.plot(bins, prof['cell_volume'])
    plt.savefig('../PigPen/test.png')


