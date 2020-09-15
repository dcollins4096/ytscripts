from go import *
import taxi
reload(taxi)
from p42_helmholtz import *
#aq20 = taxi.taxi('taxi_stops/p42_aq20.taxi')
#aq22 = taxi.taxi('taxi_stops/p42_aq22.taxi')
#aq08 = taxi.taxi('taxi_stops/p42_aq08.taxi')
#aq24 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq24_stochastic_test_1',name='p42_aq24')
#aq23 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq23_stochastic_stock',name='p42_aq23')
ca02 = taxi.load('ca02')
for t1 in [ca02]:
    for frame in [1]:
        t1.frames=[frame]
        t1.fields = ['density','x-velocity','y-velocity']
        t1.load(frame)
#t1.plot()
        field = 'acceleration'
#field = 'Driving'
        field = 'velocity'
        all_fields = ['%s-%s'%(s,field) for s in 'xyz']
#all_fields = ['DrivingField%s'%s for s in '123']
        for tmpfield in all_fields:
            fff = t1.fft(field = tmpfield)

        tests = [t1.outname]
        ntests = len(tests)
        sub_out = "_%s"*ntests%tuple(tests)
        if 0:
            MakeHelmholz(t1,frame,field)
            HelmholzPower(t1,frame,field)
            k,pd,ps=plot_helm(t1,frame,field)
            stat(pd,"PD")
            stat(ps,"PS")
        if 0:
            MakeVelocitySpectra(t1,frame)
            k, p = plot_velocity_spectra(t1,frame)
            stat(t1.ds.all_data()['y-velocity'],'VY')
        if 1:
            MakeVelocitySpectra(t1,frame)
            MakeMagneticSpectra(t1,frame)
            MakeDensitySpectra(t1,frame)
            vel_fname = spectra_filename(t1,frame,None,'velocity')
            mag_fname = spectra_filename(t1,frame,None,'magnetic')
            den_fname = spectra_filename(t1,frame,None,'density')
            mag_k,mag_p = davetools.dpy(mag_fname,['k','power'])
            vel_k,vel_p = davetools.dpy(vel_fname,['k','power'])
            den_k,den_p = davetools.dpy(den_fname,['k','power'])
            plt.close('all')
            fig,ax=plt.subplots(1,1)
            sl=slice(1,None)
            ax.plot(vel_k[sl],vel_p[sl],c='r',label='v')
            ax.plot(mag_k[sl],mag_p[sl],c='b',label='b')
            ax.plot(den_k[sl],den_p[sl],c='k',label=r'$\rho$')
            axbonk(ax,xlabel='k',ylabel='P',xscale='log',yscale='log')
            plt.savefig("%s_%04d_spectra.png"%(t1.outname,frame))
