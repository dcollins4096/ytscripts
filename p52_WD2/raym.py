
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
frame=60
rm = rainbow_map(4)
linedict={'x':'-','y':'--','z':':'}
rho_ni_label=r'$^{56}\rm{Ni}$'
flabel={'Density_56Ni':rho_ni_label,'density':r'$\rho$'}
def make_r(obj,center):
    r = np.sqrt( (obj['x']-center[0])**2 + (obj['y']-center[1])**2 + (obj['z']-center[2])**2)
    return r

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
    domain_center=  0.5*(ds.domain_left_edge+ds.domain_right_edge)
    #ray_x = ds.ortho_ray('x',[0.0,0.0])
    #argsort_x = np.argsort(ray_x['x'])
    #ray_x_field = ray_x[field][argsort_x]

    ylabel=r'$\Delta \rho $'
    ylabel=r'$\Delta \rho_{^{56}\rm{Ni}}/\bar{\rho}$'
    ylabel=r'$\Delta \rho/\bar{\rho}$'
    first_ray = None
    ax0.text(-2.0,2e8,r'$%s$'%labelmap[car.name])


    if 0:
        """rays with oray to check sanity."""
        ray_x = ds.ortho_ray('x',[0.0,0.0])
        argsort_x = np.argsort(ray_x['x'])
        ray_x_field = ray_x[field][argsort_x]
        ray_x_x = ray_x['x'][argsort_x]
        ray_r_x = make_r(ray_x, domain_center)[argsort_x]
        xok = ray_x_x>0
        ray_r_x = ray_r_x[xok]
        ray_x_x = ray_x_x[xok]
        ray_x_field = ray_x_field[xok]
        ax0.plot( ray_x_x/xscale, ray_x_field,c='k')
        for axis in 'yz':
            ray = ds.ortho_ray(axis,[0.0,0.0])
            asort = np.argsort(ray[axis])
            this_ray = ray[field][asort].v
            the_x = ray[axis][asort]
            ok_x = the_x > 0
            ok = pdf > 0
            this_r = make_r( ray, domain_center)[asort][xok]
            #this_fiducial = np.interp( the_x[ok_x], ray_x_x, ray_x_field) # this works nicely
            this_fiducial = np.interp( this_r, ray_r_x, ray_x_field) 
            delta = (this_ray[ok_x]-this_fiducial)/(0.5*(this_ray[ok_x]+this_fiducial))
            ax1.plot(the_x[ok_x]/xscale, delta, c='g')




    if 1:
        xbins, bin_center,pdf,bin_widths=toplot(profs[car.name][frame][field],quan=field)
        ok = pdf > 0
        ax0.plot( bin_center[ok]/xscale, pdf[ok],c='r')
        std = profs[car.name][frame][field].standard_deviation['gas',field]
        #ax1.plot( bin_center[ok]/xscale, std[ok]/pdf[ok]c='r')
    angles = zip( np.arange(0,np.pi/2,0.1), np.arange(0,np.pi/2,0.1))
    #angles = [[np.pi/2,0],[0,0], [np.pi/2,np.pi/2]]
    first_quan=None
    for theta,phi in angles:
        R = 2e8
        end = [R*np.sin(theta)*np.cos(phi), R*np.sin(theta)*np.sin(phi), R*np.cos(theta)]

        ray = ds.r[center:end]

        asort = np.argsort(ray['t'])
        this_quan = ray[field][asort]
        the_r = make_r(ray,domain_center)[asort]
        ax0.plot( the_r/xscale,this_quan,c=c,lw=0.1)
        if first_quan is None:
            first_quan = this_quan
            first_r = the_r
        else:
            ok = pdf > 0
            this_fiducial = np.interp( the_r, first_r, first_quan)
            print('===')
            print(ray_x_field.size)
            print(this_quan.size)

            this_delta = (this_quan.v-this_fiducial)/this_fiducial
            ax1.plot( the_r/xscale,this_delta,c=c,lw=0.1)
    #
    # oray 1: maximally good
    # 
    if 0:
        ray_x = ds.ortho_ray('x',[0.0,0.0])
        argsort_x = np.argsort(ray_x['x'])
        ray_x_field = ray_x[field][argsort_x]
        the_x = ray_x['x'][argsort_x]
        xok = the_x>0
        ax0.plot( the_x[xok]/xscale, ray_x_field[xok],c='k')
        for axis in 'yz':
            ray = ds.ortho_ray(axis,[0.0,0.0])
            asort = np.argsort(ray[axis])
            this_ray = ray[field][asort]
            the_x = ray[axis][asort]
            delta = (ray[field]-ray_x[field])/(0.5*(ray[field]+ray_x[field]))
            the_y=delta[asort]
            ax1.plot(the_x/xscale, the_y, c='k')
        
#   ylabel0=None
#   ylabel1=None
#   if col == 0:
#       ylabel0=flabel.get(field,field)
#       ylabel1=ylabel

#   axbonk(ax0,xlabel=r'',ylabel=ylabel0,yscale='log')
#   axbonk(ax1,xlabel=r'$r[10^8\ cm]$',ylabel=ylabel1,yscale='linear')
#   if ncar==0 and False:
#       yticks=ax1.get_yticks()
#       new_yticks = [expform(ytick,format="%0.1e") for ytick in yticks]
#       ax1.set_yticks(yticks)
#       ax1.set_yticklabels(new_yticks)
#   #ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

#   if not four_at_once:
#       fig.savefig('%s_%s_n%04d_radial_arms.pdf'%(field,car.name,frame))
if four_at_once:
    fig.savefig('%s_n%04d_radial_arms.pdf'%(field,frame))
    

