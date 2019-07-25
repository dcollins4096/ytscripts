from go import *
import xtra_operators as xo
reload(xo)
ip_man = "/home/dcollins/scratch/Paper52_NewDwarf/b08_uniform"
if ip_man not in sys.path:
    sys.path +=[ip_man]
    print('moo')
import initialProfile as ip
reload(ip)
car = taxi.load('p52_b08')
frame = 3
has_accel=True
def add_things(obj):
    if has_accel:
        def avec(field,data):
            output =  np.sqrt(data['Accel0']**2+data['Accel1']**2 + data['Accel2']**2)#/data['pressure'].v
            print('yes sir')
            return output
        obj.add_field('avec',avec,units='dimensionless')
        def rhoa0(field,data):
            output = -data['Accel0']*data['density'].v
            return output
        obj.add_field('rhoa0',rhoa0,units='dimensionless')
        def rhoa(field,data):
            output = data['Accel0']*data['density'].v
            return output
        obj.add_field('rhoa',rhoa,units='dimensionless')
    def dpdx(field,data):
        output = xo.grad(data,'pressure',0)
        return output
    obj.add_field('dpdx',dpdx,units='dyne/cm**3', validators=[yt.ValidateSpatial(1,['pressure'])])
    def dpdy(field,data):
        output = xo.grad(data,'pressure',1)
        return output
    obj.add_field('dpdy',dpdy,units='dyne/cm**3', validators=[yt.ValidateSpatial(1,['pressure'])])
    def dpdz(field,data):
        output = xo.grad(data,'pressure',2)
        return output
    obj.add_field('dpdz',dpdz,units='dyne/cm**3', validators=[yt.ValidateSpatial(1,['pressure'])])
    def gradmag(field,data):
        return np.sqrt( data['dpdx']**2+data['dpdy']**2+data['dpdz']**2)
    obj.add_field('gradmag',gradmag,units='dyne/cm**3', validators=[yt.ValidateSpatial(1,['pressure'])])

if 'car' not in dir() or True:
    car.operation='CenterSlice'
    car.frames=[frame]
    #car.frames=[11] #[3] #,7,11]
    car.fields=[]
    if has_accel:
        car.fields=['Accel0','Accel1','Accel2','pressure']
    car.fields+=['avec','density']#, 'pressure']#,'rhoa1','rhoa2','pressure']
    car.fields=['density','velocity_magnitude']
    car.derived_fields={'rhoa':add_things}
    car.outname = 'p52_plots/'+car.name
#   car.slice_zlim['pressure']=[8e9,4e10]
#   car.slice_zlim['avec']=[-1e12,1e12] #[1e8,1e12]
#   car.slice_zlim['avec']=[-5e10,5e10] #[1e8,1e12]
#   car.slice_zlim['Accel0']=[-5e10,5e10]
#   car.slice_zlim['Accel1']=[-5e10,5e10]
#   car.slice_zlim['Accel2']=[-5e10,5e10]
#   car.Colorbar='fixed'
#car.outname = "P52_ic/%s"%car.name

    car.axis=[0]#[0,1,2]
    #car.plot()
    ds = car.load()
    ad = ds.all_data()
if 1:
    for frame in car.frames:
        car.frames=[frame]
        ds = car.load()
        add_things(ds)
        ad = ds.all_data()
        ds.periodicity=[True,True,True]

        plt.clf()
        ext = extents()

        rho = ad['density']
        x = ad['x'].v#-0.5
        y = ad['y'].v#-0.5
        z = ad['z'].v#-0.5
        r = np.sqrt(x**2+y**2+z**2)
        grad_p_magnitude = np.sqrt( ad['dpdx']**2+ad['dpdy']**2+ad['dpdz']**2)
        if has_accel:
            ax = ad['Accel0']
            ay = ad['Accel1']
            az = ad['Accel2']
            avec = np.sqrt( ax**2+ay**2+az**2)
            plt.scatter(r,rho*avec,label=r'$\rho ||a||$',c='r',s=0.1)
            #plt.scatter(r,avec,label=r'$||a||$',c='r')
            #ext(avec.v[avec>0])
        #plt.scatter(r,rho)
        #psave(plt,'p52_plots/rho_r.png')
        plt.scatter(r,grad_p_magnitude,c='b',label='gradp',s=0.1)
        ext( grad_p_magnitude.v[ grad_p_magnitude>0])
        #plt.scatter(r,ad['pressure'],c='g',label='p')
        #plt.scatter(r, ad['TotalEnergy'],c='k',label='TE')
        #ext(grad_p_magnitude.v)
        #ext(ad['pressure'].v)
        #ext(ad['TotalEnergy'].v)
        #plt.ylim(ext.minmax)
        #plt.yscale('symlog',linthreshy = 1)
        if 0:
            rr=ip.r[0]+0
            pp=ip.P[0]+0
            rmid = 0.5*(rr[1:]+rr[:-1])[1:]
            dr = (rr[1:]-rr[:-1])[1:]
            dp  = (pp[1:]-pp[:-1])[1:]
            dpdr=np.abs(dp/dr)
            theg=np.abs(ip.g)*ip.rho[0]
            plt.plot( rmid, dpdr, c='b')
            plt.plot( rr, theg, c='g')
        plt.yscale('symlog',simlogthreshy=1e1)
        plt.ylim( ext.minmax)
        #plt.yscale('log')
        plt.legend(loc=0)
        psave(plt,'p52_plots/dpdx_r_%04d.png'%frame)

        if 1:
            fig,ax=plt.subplots(1,1)
            ax.scatter(r, rho,label='rho')
            axbonk(ax,xlabel='r',ylabel='rho',xscale='log',yscale='log')
            psave(fig,'p52_plots/rho_scatter_%04d.png'%frame)
            ax.clear()
            ax.scatter(r,ad['velocity_magnitude'],label='||v||')
            psave(fig,'p52_plots/vel_scatter_%04d.png'%frame)
            plt.close(fig)
            


        plt.clf()
        ex = extents()
        p = np.abs(ad['pressure']+0)
        #p[ r==r.min()]=p.min()
        plt.scatter(r,p,c='b',label='p')
        #if has_accel:
        #    plt.scatter(r,avec,label=r'$||a||$',c='r')
        #    ex(avec.v)

        plt.xlabel(r'$r$')
        plt.ylabel(r'$acceleration$')
        ex(np.abs(p.v))

        if 0:
            plt.plot( rr, pp, c='k')


        plt.ylim(ex.minmax)
        plt.legend(loc=0)
        #plt.plot(r,350*r)
        plt.yscale('log')
        plt.plot(r,r)
        psave(plt,'p52_plots/pressure_%04d.png'%frame)


