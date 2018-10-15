if 'ef' not in dir():
    execfile('go')
import taxi
reload(taxi)

execfile('p33_sims.py')
all_fields = all_fields_1 + ['Metal_Density']
#all_fields = all_fields_simple
L=len(all_fields)
phase_list = zip(['density']*L,all_fields)
#phase_list += [['density','magnetic_field_strength']]
#phase_list += [['Temperature','magnetic_field_strength']]
fields = np.unique(nar(phase_list).flatten())

if 'aj01' not in dir():
    aj01=taxi.taxi('aj01.taxi')
    aj02=taxi.taxi('aj02.taxi')
    aj07=taxi.taxi('aj07.taxi')
    aj08=taxi.taxi('aj08.taxi')
    f1 = taxi.fleet([aj07,aj08])
    taxi_list = [aj01, aj02] #, aj05]
    execfile('p33_sims.py')
    all_fields = all_fields_1 + ['Metal_Density']
    #all_fields = all_fields_simple
    L=len(all_fields)
    phase_list = zip(['density']*L,all_fields)
    #phase_list += [['density','magnetic_field_strength']]
    #phase_list += [['Temperature','magnetic_field_strength']]
    fields = np.unique(nar(phase_list).flatten())
    if 0:
        for car in taxi_list:
            car.fields = ['density','Temperature']
            car.frames=range(0,100,2)
            car.frames = [30,40]
            #val,car.center=car.ds.find_max('density')
            car.center= [0.5]*3 #car.ds.arr([0.5]*3,'code_length')
            car.radius = 0.1
            car.operation='RegionProjection'
            car.region_type='sphere'
            car.width=(600,'kpc')
            car.restrict=True
            #car.plot()

if 0:
    reload(taxi)
    aj01=taxi.taxi('aj01.taxi')
    aj02=taxi.taxi('aj02.taxi')
    taxi_list = [aj01,aj02]
    for car in taxi_list:

        car.operation='RegionProjection'
        car.region_type='sphere'
        car.callbacks=['particles']
        car.width=(600,'kpc')
        car.restrict=True
        car.frames = range(40,51)
        car.axis='x'
        car.plot()

if 0:
    for n,car in enumerate(taxi_list):
        fields = ['TotalEnergy']
        car.frames = [190]
        if n>0:
            car.extrema=last_extrema
        car.find_extrema(fields)
        last_extrema=car.extrema
    for n,car in enumerate(taxi_list):
        car.extrema=last_extrema

if 0:
    L=len(all_fields)
    phase_list = phase_list = zip(['density']*L,all_fields)
    for n,car in enumerate(taxi_list):
        car.Colorbar = 'fixed'
        for x,y in phase_list:
            car.phase([x,y,'cell_mass'])

if 0:
    index_list = ['aj01','aj02']
    particles={}
    particles['aj01']={}
    car = aj01
    for n in range(100):
        print "N %04d"%n,
        car.fill(n)
        if car.ds['NumberOfParticles'] > 0:
            iii = [int(mf) for mf in sorted(car.ds.all_data()['particle_index'].v)]
            particles[car.name][n]=iii
            last_p = int(iii[0])-1
            for p in iii:
                p = int(p)
                difference = p - last_p - 0
                for n in range(difference):
                    print "  ",
                print p,
                last_p = p
        else:
            print "--",
        print ""


if 0:
    for car in taxi_list:
        car.extrema = {'Temperature': [np.array(2137.9392148662896), np.array(136413446.8401202)], 'density': [np.array(4.8094553109589204e-33), np.array(3.413387818580601e-24)]}
        car.Colorbar = 'fixed'
        car.phase(['density','Temperature','cell_mass'])

if 0:
    for car in taxi_list:
        car.fields = ['density']
        car.Colorbar = 'monotonic'
        car.operation='RegionProjection'
        car.region_type='sphere'
        car.axis=['y']
        car.restrict=True
        car.plot()
if 0:
    frame=600
    fname = '/mnt/c/scratch/sciteam/dcollins/Paper33_Galaxy/aj01_fiducial/DD%04d/DD%04d'%(frame,frame)
    aj01.fill(0)
    ds = aj01.ds
    sphere = ds.sphere([0.5]*3,0.1)
    proj = ds.proj('Temperature',0, data_source = sphere)
    pw=proj.to_pw()
    pw.set_width(600,'kpc')
    print pw.save('aj01_t7')

if 0:
    frame=600
    fname = '/mnt/c/scratch/sciteam/dcollins/Paper33_Galaxy/aj01_fiducial/DD%04d/DD%04d'%(frame,frame)
    ds = yt.load(fname)
    sphere = ds.sphere
    proj = ds.proj('Temperature',0)
    pw=proj.to_pw()
    pw.set_width(600,'kpc')
    print pw.save('aj01_t1')
