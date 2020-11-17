from go import *

car = taxi.load('p52_434')
car.frames=[0]
field='magnetic_field_strength'
if 0:
    car.profile([field,'cell_volume'])
if 'ds' not in dir():
    ds = car.load()
    ad=ds.all_data()
    nbins=64
if 'prof' not in dir():
    prof=yt.create_profile(ad,[field],fields=['cell_volume'],weight_field=None,n_bins=nbins,fractional=True) #override_bins=override_bins)

if 'radial' not in dir():
    radial=yt.create_profile(ad,['radius'],fields=[field],weight_field='cell_volume',n_bins=nbins) #override_bins=override_bins)
if 0:
    bins = prof.x_bins
    bincen = 0.5*(bins[:-1]+bins[1:])
    fig,ax=plt.subplots(1,1)
    ax.plot(bincen, prof['cell_volume'])
    axbonk(ax,xlabel='B',ylabel='P(B)db',xscale='log',yscale='log')
    fig.savefig('pdf_test.png')
if 1:
    rads = radial.x_bins
    bincen = 0.5*(rads[:-1]+rads[1:])
    fig,ax=plt.subplots(1,1)
    ax.plot(bincen,radial[field])
    axbonk(ax,xlabel='r[cm]', ylabel=field,xscale='log',yscale='log')
    fig.savefig('bvsr.png')
    plt.close(fig)


plt.close(fig)
