if 'ef' not in dir():
    execfile('go')
import taxi
reload(taxi)
#t1 = taxi.taxi('taxi_stops/p42_aq20.taxi')
#t1 = taxi.taxi('taxi_stops/p42_aq22.taxi')
t1 = taxi.taxi('taxi_stops/p42_aq08.taxi')
t1 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq24_stochastic_test_1',name='p42_aq24')
frame = 40
t1.frames=[frame]
t1.fields = ['density']
#t1.plot()
field = 'acceleration'
#field = 'Driving'
#field = 'velocity'
all_fields = ['%s-%s'%(s,field) for s in 'xyz']
#all_fields = ['DrivingField%s'%s for s in '123']
for tmpfield in all_fields:
    fff = t1.fft(field = tmpfield)

tests = [t1.outname]
ntests = len(tests)
sub_out = "_%s"*ntests%tuple(tests)
ef('p42_helmholtz.py')
MakeHelmholz(t1,frame,field)
HelmholzPower(t1,frame,field)
k,pd,ps=plot_helm(t1,frame,field)
