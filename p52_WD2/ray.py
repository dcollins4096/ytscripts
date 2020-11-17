
from go import *
reload(taxi)
#import p52_WD2.multiplot2 as mp2
#reload(mp2)
cars = ['p52_441','p52_432','p52_433','p52_434']
labelmap = {'p52_441':'D22','p52_432':'D26','p52_433':'D27','p52_434':'D28'}
if 'flt' not in dir():
    flt = taxi.fleet(cars)

four_at_once = True
if four_at_once:
    fig, ax=plt.subplots(4,2,sharex=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    #fig, axlist, cb = mp2.mp(
    # GOT HERE
#field='Density_56Ni'
#field='magnetic_field_strength'
field='density'
frame=50
rm = rainbow_map(4)
linedict={'x':'-','y':'--','z':':'}
rho_ni_label=r'$^{56}\rm{Ni}$'
flabel={'Density_56Ni':rho_ni_label,'density':r'$\rho$'}
for ncar,car in enumerate(flt.taxi_list):
    c=rm(ncar)
    print(c)
    if four_at_once:
        row = ncar//2
        col = ncar%2
        ax0 = ax[2*row  ][col]
        ax1 = ax[2*row+1][col]
        if col:
            ax0.yaxis.tick_right()
            ax1.yaxis.tick_right()

        special_y_axis=False
        if ncar==0:
            special_y_axis=True
        xscale = 1e8
    else:
        ax0.clear()
        ax1.clear()
    ds = car.load(frame)
    ray_x = ds.ortho_ray('x',[0.0,0.0])
    argsort_x = np.argsort(ray_x['x'])
    ray_x_field = ray_x[field][argsort_x]

    ylabel=r'$\Delta \rho $'
    ylabel=r'$\Delta \rho_{^{56}\rm{Ni}}/\bar{\rho}$'
    ylabel=r'$\Delta \rho/\bar{\rho}$'
    for axis in 'yz':
        ray = ds.ortho_ray(axis,[0.0,0.0])
        asort = np.argsort(ray[axis])
        this_ray = ray[field][argsort]
        the_x = ray[axis][asort]/xscale
        delta = (ray[field]-ray_x[field])/(0.5*(ray[field]+ray_x[field]))
        the_y=delta[asort]
        if special_y_axis:
            the_y /= 1e-4
            ax1.text(-2,-5,r'$\times 10^{-4}$')
        ax0.plot( the_x,ray[field][asort],c=c,linestyle=linedict[axis])
        ax1.plot( the_x,the_y,c=c,linestyle=linedict[axis])
        ax0.text(-2.0,2e8,r'$%s$'%labelmap[car.name])
    asort=np.argsort(ray_x['x'])
    ylabel0=None
    ylabel1=None
    if col == 0:
        ylabel0=flabel.get(field,field)
        ylabel1=ylabel

    ax0.plot( the_x, ray_x_field,c=c,linestyle=linedict['x'])
    axbonk(ax0,xlabel=r'',ylabel=ylabel0,yscale='log')
    axbonk(ax1,xlabel=r'$r[10^8\ cm]$',ylabel=ylabel1,yscale='linear')
    if ncar==0 and False:
        yticks=ax1.get_yticks()
        new_yticks = [expform(ytick,format="%0.1e") for ytick in yticks]
        ax1.set_yticks(yticks)
        ax1.set_yticklabels(new_yticks)
    #ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

    if not four_at_once:
        fig.savefig('%s_%s_n%04d_radial_arms.pdf'%(field,car.name,frame))
if four_at_once:
    fig.savefig('%s_n%04d_radial_arms.pdf'%(field,frame))
    

