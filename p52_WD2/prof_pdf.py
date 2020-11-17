
from go import *
reload(taxi)
cars = ['p52_441','p52_432','p52_433','p52_434']
if 'flt' not in dir():
    flt = taxi.fleet(cars)
flt=taxi.fleet(['p52_434'])
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

    use_delta = False
    rm = rainbow_map(4)
    frame_list=[0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100] 
    frame_list = [0,10,50]
    frame_list=[50]
    frame_list=[0]
    line_list={0:'-',60:'--'}
    line_list={0:'-',1:'--',2:':',3:'-.'}
    line_list={0:'-',10:'--',50:'-.'}
    field_list = ['magnetic_field_strength']; prefix='mag_pdf'
    leg=dict(zip(['Density_56Ni','magnetic_energy','thermal_energy','kinetic_energy','pressure','density','velocity_magnitude','magnetic_field_strength'],['Ni','BE','TE','KE','P','rho','v','B']))
    

    if 'pdfs' not in dir():
        pdfs = {}
    if 1:
        for car in flt.taxi_list:
            if car.name not in pdfs:
                pdfs[car.name]={}
            for frame in frame_list:
                if frame not in pdfs[car.name]:
                    pdfs[car.name][frame]={}
	
        ds = car.load()
        nbins=64
        #rmax = (3/4)**0.5*ds.domain_width[0] #half the diagonal
        #rmin = rmax**(1/nbins)
        #override_bins={'radius':np.logspace(np.log10(rmin),np.log10(rmax),16)}


        for car in flt.taxi_list: 
            for frame in frame_list:
                ds=car.load(frame)
                ad=ds.all_data()
                for field in field_list:
                    if field not in pdfs[car.name][frame]:
                        prof=yt.create_profile(ad,[field],fields=['cell_volume'],weight_field=None,n_bins=nbins,fractional=True) #override_bins=override_bins)
                        pdfs[car.name][frame][field]=prof

    if use_delta:
        fig, axes = plt.subplots(2,1,sharex=True)
        fig.subplots_adjust(wspace=0, hspace=0)
        ax=axes[0]
        ax1=axes[1]
    else:
        fig,ax = plt.subplots(1)
    for nframe, frame in enumerate(frame_list):
        ax.clear()
        for ncar,car in enumerate(flt.taxi_list): 
            for nf,field in enumerate(field_list):

                xbins, bin_center,pdf,bin_widths=toplot(pdfs[car.name][frame][field],quan='cell_volume')
                if field in ['pressure']:
                    c='k'
                else:
                    c=rm(ncar)
                label = '%s %s'%(car.name[4:],leg[field])
                ax.plot(bin_center,pdf,c=c,linestyle=line_list.get(nf,':'), label=label)
                ax.set_xscale('log');ax.set_yscale('log')
                fig.savefig('temp.png')

                if use_delta:
                    xbins0, bin_center0,pdf0,bin_widths0=toplot(pdfs['p52_441'][frame][field],quan='cell_volume')
                    delta = (pdf - pdf0)/(0.5*(pdf+pdf0))
                    ax1.plot(bin_center,delta,c=c,linestyle=line_list.get(nf,':'), label=label)
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
        elif prefix == 'mag_pdf':
            yscale= [1e-6,2e-1]
            yscale=[prof['cell_volume'].min(), prof['cell_volume'].max()]
        else:
            yscale = ext.minmax
        this_map={'magnetic_field_strength':'B'}
        ylabel='P(%s)'%this_map.get(field,field)
        xlabel=this_map.get(field,field)
        axbonk(ax,xscale='log',yscale='log',ylim=yscale, ylabel=ylabel, xlabel=xlabel)
        #ax.legend(loc=0)
        fig.savefig('%s_n%04d.pdf'%(prefix,frame))
        plt.close('all')

