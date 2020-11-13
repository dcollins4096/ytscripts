from go import *
reload(taxi)
if 'flt' not in dir():
    flt = taxi.fleet(['p52_441','p52_432','p52_433','p52_434'])
if 0:
    flt['fields'] = ['density','Density_56Ni']
    flt['frames'] = [60]
    flt['callbacks'] = ['magnetic_streamlines']
    flt.profile(['radius','Density_56Ni'], weight_field='cell_volume')
if 0:
    flt['callback_args']={'time_title':{'units':'s','format':'%0.4f'}}
    flt['callbacks']=['time_title']
    flt['frames']=[60]
    flt['fields']=['Density_56Ni']
    flt.plot()
if 0:
    flt['frames'] = [60]
    flt['fields'] = ['pressure']
    flt.plot()
if 0:
    flt['frames']=[0]
    flt.profile(['radius','TotalEnergy'],weight_field='cell_volume')

if 0:
    flt['frames']=[0]
    flt.profile(['radius','magnetic_energy'],weight_field='cell_volume')

if 0:
    fig,ax = plt.subplots(1,1)
    ex = extents()
    for car in flt.taxi_list:
        bins = car.profile_data['all_xbins'][0]
        proo = car.profile_data['all_profiles'][0]
        print('wuut',proo.min(),proo.max())
        ex(proo)
        ax.plot(bins,proo,label=car.name)
    axbonk(ax,xlabel=car.profile_data['fields'][0],ylabel=car.profile_data['fields'][1],xscale='log',yscale='log',ylim=ex.minmax)
    fig.savefig('thing2.png')
    plt.close('all')

if 1:
    for car in flt.taxi_list:
        car.profile(['radius','magnetic_field_strength'],weight_field='cell_volume')

if 0:
    fig,ax = plt.subplots(1,1)
    ex = extents()
    for car in flt.taxi_list:
        ds=car.load(0)
        ad=ds.all_data()
        prof=yt.create_profile(ad,['radius'],fields=['magnetic_field_strength'])
        xbinsc =0.5*( prof.x_bins[1:] + prof.x_bins[:-1])
    axbonk(ax,xscale='log',yscale='log')
    fig.savefig('thing3.png')
    plt.close('all')



if 0:
    fig, ax = plt.subplots(1,1)
    for n,car in enumerate(flt.taxi_list):
        if n == 2:
            continue
        the_x = car.profile_data['all_xbins'][0] * (1+0.1*n)

        the_y = car.profile_data['all_profiles'][0]
        labels = car.profile_data['fields']
        ax.plot(the_x,the_y,label=car.name)
    axbonk(ax,xlabel=labels[0],ylabel=labels[1],yscale='log',xscale='log')
    ax.legend(loc=0)
    outname='plots_to_sort/total_energy.pdf'
    fig.savefig(outname)
    print(outname)


    plt.close(fig)

#flt.profile(['radius','density'], weight_field='cell_volume')
