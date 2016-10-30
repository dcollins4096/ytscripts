if 'ef' not in dir():
    execfile('go')
import taxi
reload(taxi)
aq20 = taxi.taxi('taxi_stops/p42_aq20.taxi')
aq22 = taxi.taxi('taxi_stops/p42_aq22.taxi')
aq08 = taxi.taxi('taxi_stops/p42_aq08.taxi')
aq24 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq24_stochastic_test_1',name='p42_aq24')
aq23 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq23_stochastic_stock',name='p42_aq23')
for t1 in [aq22]:
    for frame in [0]:
        t1.frames=[frame]
        t1.fields = ['density','x-velocity','y-velocity']
        t1.fill(frame)
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
        ef('p42_helmholtz.py')
        if 1:
            MakeHelmholz(t1,frame,field)
            HelmholzPower(t1,frame,field)
            k,pd,ps=plot_helm(t1,frame,field)
            stat(pd,"PD")
            stat(ps,"PS")
        if 0:
            MakeVelocitySpectra(t1,frame)
            k, p = plot_velocity_spectra(t1,frame)
            stat(t1.ds.all_data()['y-velocity'],'VY')
