
from go import *
reload(taxi)
cars = ['p52_441','p52_432','p52_433','p52_434']
if 'flt' not in dir():
    flt = taxi.fleet(cars)
#car_list=[]
#for name in cars:
#    car_list.append(taxi.load(name))

def toplot(prof,quan = 'cell_volume'):
    xbins = prof.x_bins
    bin_center = 0.5*(xbins[1:]+xbins[:-1])
    bin_widths = xbins[1:]-xbins[:-1]
    pdf = prof[quan]
    pdf = pdf/bin_widths
    return xbins, bin_center,pdf,bin_widths

if 'ext' not in dir():
    ext = extents()
if 1:

    rm = rainbow_map(4)
    frame_list=[0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100] 
    frame_list = [0,1,2,3,10]
    line_list={0:'-',60:'--'}
    line_list={0:'-',1:'--',2:':',3:'-.'}
    #field_list=['Density_56Ni']; prefix='nickle_profiles'
    #field_list=['magnetic_energy','kinetic_energy','pressure']; prefix='energy_profiles'
    #field_list=['magnetic_energy','thermal_energy','kinetic_energy','pressure']; prefix='energy_profiles'
    field_list = ['density']; prefix='density_profiles'
    #field_list = ['magnetic_energy','pressure']; prefix='magnetic_profiles'
    #field_list=['velocity_magnitude']; prefix='velocity_profile'
    leg=dict(zip(['Density_56Ni','magnetic_energy','thermal_energy','kinetic_energy','pressure','density','velocity_magnitude'],['Ni','BE','TE','KE','P','rho','v']))
    

    if 'profs' not in dir():
        profs = {}
    if 1:
        for car in flt.taxi_list:
            if car.name not in profs:
                profs[car.name]={}
            for frame in frame_list:
                if frame not in profs[car.name]:
                    profs[car.name][frame]={}
	
        ds = car.load()
        nbins=16
        #rmax = (3/4)**0.5*ds.domain_width[0] #half the diagonal
        #rmin = rmax**(1/nbins)
        #override_bins={'radius':np.logspace(np.log10(rmin),np.log10(rmax),16)}


        for car in flt.taxi_list: 
            for frame in frame_list:
                ds=car.load(frame)
                ad=ds.all_data()
                for field in field_list:
                    if field not in profs[car.name][frame]:
                        prof=yt.create_profile(ad,['radius'],fields=[field],n_bins=nbins) #override_bins=override_bins)
                        profs[car.name][frame][field]=prof

    fig,ax = plt.subplots(1)
    for nframe, frame in enumerate(frame_list):
        ax.clear()
        for ncar,car in enumerate(flt.taxi_list): 
            for nf,field in enumerate(field_list):

                xbins, bin_center,pdf,bin_widths=toplot(profs[car.name][frame][field],quan=field)
                if field in ['pressure']:
                    c='k'
                else:
                    c=rm(ncar)
                label = '%s %s'%(car.name[4:],leg[field])
                ax.plot(bin_center,pdf,c=c,linestyle=line_list.get(nf,':'), label=label)
                if np.abs(pdf).sum() > 0:
                    ext(pdf[pdf>0].v)
                print(nframe)


        if 1:
            time = 1.*frame/100.
            print("TIIIIIIIIM", time)
            ax.set_title(r"$t=%0.2f\ \rm{s}$"%time)
        ex = extents()
        yscale=None
        if prefix == 'velocity_profile':
            yscale=[1e-5,3e2]
        elif prefix == 'energy_profiles':
            #for energy
            yscale=[5e2,5e24]
        else:
            yscale = ext.minmax
        axbonk(ax,xscale='log',yscale='log',ylim=yscale, ylabel='Energy', xlabel='r')
        #ax.legend(loc=0)
        fig.savefig('%s_n%04d.png'%(prefix,frame))
        plt.close('all')

