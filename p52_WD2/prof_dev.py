
from go import *
reload(taxi)
cars = ['p52_441']#,'p52_432','p52_433','p52_434']
if 'flt' not in dir():
    flt = taxi.fleet(cars)


def toplot(prof,quan = 'cell_volume'):
    xbins = prof.x_bins
    bin_center = 0.5*(xbins[1:]+xbins[:-1])
    bin_widths = xbins[1:]-xbins[:-1]
    pdf = prof[quan]
    pdf = pdf/bin_widths
    return xbins, bin_center,pdf,bin_widths


frame_list=[0]

field_list=['Density_56Ni'];prefix='Ni'
field_list=['density'];prefix='rho'
#field_list=['pressure'];prefix='Pressure'

def add_squares(obj):
    def square_ni(field,data):
        return data['Density_56Ni']**2
    obj.add_field('Density_56Ni_squared',function=square_ni)
    def square_p(field,data):
        return data['pressure']**2
    obj.add_field('pressure_squared',function=square_p, units='dyne**2/cm**4',sampling_type='cell')
    def square_rho(field,data):
        return data['density']**2
    obj.add_field('density_squared',function=square_rho, units='g**2/cm**6',sampling_type='cell')

if 'profs' not in dir():
    profs={}
    profs2={}
    nbins=16

weight_field='cell_volume'
if 1:
    for car in flt.taxi_list:
        if car.name not in profs:
            profs[car.name]={}
            profs2[car.name]={}
        for frame in frame_list:
            if frame not in profs[car.name]:
                profs[car.name][frame]={}
                profs2[car.name][frame]={}

    for car in flt.taxi_list: 
        for frame in frame_list:
            ds=car.load(frame)
            add_squares(ds)
            ad=ds.all_data()
            for field in field_list:
                if field not in profs[car.name][frame]:
                    prof2=yt.create_profile(ad,['radius'],fields=[field+"_squared"],n_bins=nbins,weight_field=weight_field) #override_bins=override_bins)
                    profs2[car.name][frame][field]=prof2

                    prof=yt.create_profile(ad,['radius'],fields=[field],n_bins=nbins,weight_field=weight_field) #override_bins=override_bins)
                    profs[car.name][frame][field]=prof

if 1:
    nbin = 0
    mask=(ad['radius'] >= profs[car.name][frame][field].x_bins[nbin])*(ad['radius'] < profs[car.name][frame][field].x_bins[nbin+1])
    rho_bin = ad[field][mask]
    dv_bin = ad['cell_volume'][mask].sum()
    my_mean = (rho_bin*dv_bin).sum()/dv_bin.sum()
    xbins, bin_center,pdf,bin_widths=toplot(profs[car.name][frame][field],quan=field)

    print("%0.3e %0.3e %0.3e"%(pdf[nbin-1],pdf[nbin],pdf[nbin+1]))
    print("     %0.3e"%my_mean)


if 1:
    fig,ax = plt.subplots(1)
    for nframe, frame in enumerate(frame_list):
        ax.clear()
        for ncar,car in enumerate(flt.taxi_list): 
            for nf,field in enumerate(field_list):

                xbins, bin_center,pdf,bin_widths=toplot(profs[car.name][frame][field],quan=field)
                xbins, bin_center,pdf2,bin_widths=toplot(profs2[car.name][frame][field],quan=field+"_squared")
                #if field in ['pressure']:
                #    c='k'
                #else:
                #    c=rm(ncar)
                c=None
                label = car.name #'%s %s'%(car.name[4:],leg[field])
                line_list={}
                #ax.errorbar(bin_center,pdf.v,c=c,linestyle=line_list.get(nf,':'), label=label,yerr=np.sqrt(pdf2.v-pdf.v**2))
                ax.plot(bin_center,pdf.v**2,c=c,linestyle=line_list.get(nf,':'), label=label)
                ax.plot(bin_center,pdf2.v,c=c,linestyle=line_list.get(nf,':'), label=label)
                #if np.abs(pdf).sum() > 0:
                #    ext(pdf[pdf>0].v)
                #print(nframe)

    axbonk(ax,xlabel='r',ylabel=field,xscale='log',yscale='log')
    ax.legend(loc=1)
    fig.savefig('deviation_%s_n%04d.pdf'%(field,frame))
    plt.close(fig)
