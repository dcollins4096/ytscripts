
from go import *
reload(taxi)
cars = ['p52_441','p52_432','p52_433','p52_434']
if 'flt' not in dir():
    flt = taxi.fleet(cars)
labelmap = {'p52_441':'D22','p52_432':'D26','p52_433':'D27','p52_434':'D28'}
#car_list=[]
#for name in cars:
#    car_list.append(taxi.load(name))

def toplot(prof,quan = 'cell_volume'):
    xbins = prof.x_bins
    bin_center = 0.5*(xbins[1:]+xbins[:-1])
    bin_widths = xbins[1:]-xbins[:-1]
    pdf = prof[quan]
    pdf = pdf#/bin_widths
    return xbins, bin_center,pdf,bin_widths

if 'ext' not in dir():
    ext = extents()
if 1:

    use_delta = False
    four_at_once = False
    y_tick_right=True
    rm = rainbow_map(4)
    frame_list=[0]
    line_list={0:'-',60:'--'}
    line_list={0:'-',1:'--',2:':',3:'-.'}
    line_list={0:'-',1:'-'}
    field_list = ['pressure','magnetic_energy']; prefix='pressure_profiles'
    leg=dict(zip(['Density_56Ni','magnetic_energy','thermal_energy','kinetic_energy','pressure','density','velocity_magnitude','magnetic_field_strength'],['Ni','BE','TE','KE','P','rho','v','B']))
    

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
        for car in flt.taxi_list: 
            for frame in frame_list:
                ds=car.load(frame)
                ad=ds.all_data()
                for field in field_list:
                    if field not in profs[car.name][frame]:
                        prof=yt.create_profile(ad,['radius'],fields=[field],n_bins=nbins) #override_bins=override_bins)
                        profs[car.name][frame][field]=prof

    if use_delta:
        fig, axes = plt.subplots(2,1,sharex=True, figsize=(2,2))
        fig.subplots_adjust(wspace=0, hspace=0)
        ax=axes[0]
        ax1=axes[1]
    elif four_at_once:
        fig, axes = plt.subplots(2,2,sharex=True,sharey=True)
        fig.subplots_adjust(wspace=0, hspace=0)
        ax_list=axes.flatten()
    else:
        #fig,ax = plt.subplots(1,figsize=(4,3))
        fig,ax = plt.subplots(1)
    if y_tick_right:
        ax_twin = ax.twinx()

    ex = extents()
    for nframe, frame in enumerate(frame_list):
        if  four_at_once:
            ax = ax_list[nframe]
        else:
            ax.clear()
        for ncar,car in enumerate(flt.taxi_list): 
            for nf,field in enumerate(field_list):

                xbins, bin_center,pdf,bin_widths=toplot(profs[car.name][frame][field],quan=field)
                if field in ['pressure']:
                    c='k'
                    label=None
                    if car.name in ['p52_441']:
                        label=r'$P_{\rm{gas}} \rm{(all)}$'
                else:
                    c=rm(ncar)
                    label=r'$B^2\ \rm{%s}$'%labelmap[car.name] 
                linestyle=line_list.get(nf,':')
                ax.plot(bin_center,pdf,c=c,linestyle=linestyle, label=label)
                if pdf.sum() > 0:
                    ex(pdf[pdf>1e-16].v)
                if use_delta:
                    xbins0, bin_center0,pdf0,bin_widths0=toplot(profs['p52_441'][frame][field],quan=field)
                    delta = (pdf - pdf0)/(0.5*(pdf+pdf0))
                    ax1.plot(bin_center,delta,c=c,linestyle=line_list.get(nf,':'), label=label)
                if np.abs(pdf).sum() > 0:
                    ext(pdf[pdf>0].v)
                print(nframe)


        if not four_at_once:
            time = 1.*frame/100.
            print("TIIIIIIIIM", time)
            ax.set_title(r"$t=%0.2f\ \rm{s}$"%time)
        else:
            time = 1.*frame/100.
            ax.text(5e6,1e10,r"$t=%0.2f\ \rm{s}$"%time)
        yscale=None
        ylabel='PDF'
        xlabel=r'$r\ \rm{[cm]}$'
        if prefix == 'velocity_profile':
            ylim=[1e-5,3e2]
        elif prefix == 'energy_profiles':
            pass
            #for energy
            #ylim=[5e2,5e24]
            ylim=ext.minmax
            ylim=[1e6,5e29]
            if nframe%2==0:
                ylabel=r'$E\ [\rm{g}\rm{cm}^{-3}]$'
            else:
                ylabel=""
            print(ylabel)
        elif prefix == 'mag_profile':
            ylim= None
            ylabel= r'$B\ \rm{[G]}$'
        elif prefix == 'pressure_profiles':
            ylim= [1e9,1e30]
            ylabel=r'$\rm{Pressure} [\rm{erg}/\rm{cm^3}]$'
            ax.legend(loc=0)
        else:
            ylim = ext.minmax

        axbonk(ax,xscale='log',yscale='log',ylim=ylim, ylabel=ylabel, xlabel=xlabel)
        ax.set_ylim(ylim)
        if y_tick_right:
            ylim = ax.get_ylim()
            ax_yticks=ax.get_yticks()
            new_yticks= ax_yticks**0.5
            ax_twin.set_yticks(new_yticks)
            y0,y1=ax.get_ylim()
            ax_twin.set_ylim( (8*np.pi*y0)**0.5, (8*np.pi*y1)**0.5)
            ax_twin.set_yscale('log')
            ax_twin.set_ylabel(r'$B\ \rm{[G]}$')
        fig.savefig('%s_n%04d.pdf'%(prefix,frame))
        plt.close('all')

