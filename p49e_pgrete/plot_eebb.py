
from go import *
import turb_quan
reload(turb_quan)
flt_list=['pg%d'%s for s in range(12)]
flt = taxi.fleet(flt_list)

for car in flt.taxi_list:
    qb = turb_quan.quan_box(car)
    for frame in range(10,110,10):
        qb.EBall([frame])
        plt.clf()
        for nax,ax in enumerate('xyz'):
            clfile = "%s/frbs/DD%04d_Cl%s.dat"%(car.directory,frame,ax)
            b = np.loadtxt(clfile)
            ell = b[:,0]
            EE = b[:,1]
            BB = b[:,2]
            ok = BB>0
            c='rgb'[nax]
            do=1
            if do==0:
                plt.plot( ell[ok], EE[ok]/BB[ok],c=c)
                outname = "%s_DD%04d_ratio_take1.png"%(car.outname,frame)
                plt.ylim(1,5)
                plt.yscale('linear')
            if do==1:
                plt.plot( ell, EE*ell**3.8, c=c,linestyle='-')
                plt.plot( ell, BB*ell**3.8, c=c,linestyle=':')
                outname = "%s_DD%04d_compensated_take1.png"%(car.outname,frame)
                ylabel = r'$T_{xx} \ell^{3.8}$'
                plt.ylim(5e2,1e4)
                plt.yscale('log')
        print(outname)
        plt.xlabel(r'$\ell$')
        plt.xscale('log')
        plt.ylabel(ylabel)
        plt.savefig(outname)

        #qb.QUEB(frame)
        #qb.EBslopes(frame)
