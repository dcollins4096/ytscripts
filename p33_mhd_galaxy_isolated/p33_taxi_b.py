if 'ef' not in dir():
    execfile('go')
import taxi
reload(taxi)

aj01=taxi.taxi(directory='/Users/dcollins/scratch/Paper33_galaxies/aj01',name='aj01')
aj02=taxi.taxi(directory='/Users/dcollins/scratch/Paper33_galaxies/aj02',name='aj02')
taxi_list = [aj01, aj02]
for car in taxi_list:
    car.fields = ['density','Temperature']
    car.frames=range(30,40)
    car.fill(0)
    #val,car.center=car.ds.find_max('density')
    car.center= [0.5]*3 #car.ds.arr([0.5]*3,'code_length')
    car.radius = 0.1
    car.operation='RegionProjection'
    car.region_type='sphere'
    car.width=(600,'kpc')
    car.axis=0
    car.restrict=True
    car.plot()
    

if 0:
    for car in taxi_list:
        car.find_extrema()
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
