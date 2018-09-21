
execfile('go')
execfile('p42_pg.py')
execfile('p42_spectra.py')
ef('p42_quan.py')
import taxi

all_taxi=['aw05', 'aw06', 'aw07', 'aw08', 'aw09', 'aw10',  'aw15', 'aw16', 'aw17', 'aw18']
fields = ['%s-velocity'%s for s in 'xyz']
for carname in all_taxi:
    car = taxi.taxi(carnmae)
    car.fill(0)
    all_frames = car.frame_dict.keys()
    car.frames = all_frames[10::10]
    if all_frames[-1] not in car.frames:
        car.frames += [all_frames[-1]]
    for frame in car.frames:
        outname = 'p42_power_%s_%04d.pdf'%(car.name, frame)
        if len(glob.glob(outname) ) > 0:
            continue
        all_the_spectra(car,fields)
        sys.exit(0)
