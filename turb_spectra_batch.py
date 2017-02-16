"""
Spectra generator.
python turb_spectra_batch.py taxi_name frames
Frames is  optional.
If frames='car', all the frames on the taxi file are used.
If frames is a list, the list is used.
Otherwise, every tenth frame plus the last frame is used.

Skips all extant spectra.
Exits upon completion of one spectra for memory issues (python gc is not good for extreme memory usage.)
"""
execfile('go')
execfile('p42_pg.py')
execfile('p42_spectra.py')
import taxi

carname=sys.argv[1] #['aw05', 'aw06', 'aw07', 'aw08', 'aw09', 'aw10',  'aw15', 'aw16', 'aw17', 'aw18']
fields = ['%s-velocity'%s for s in 'xyz']
car = taxi.taxi(carname)
frames=None
if len(sys.argv) > 2:
    if sys.argv[2] == 'car':
        frames = car.frames
    else:
        frames = [int(a) for a in sys.argv[2:]]

if frames is None:
    all_frames = car.frame_dict.keys()
    car.frames = all_frames[10::10]
    if all_frames[-1] not in car.frames:
        car.frames += [all_frames[-1]]
else:
    car.frames=frames
car.frames = car.frames[::-1]
for frame in car.frames:
    outname = 'p42_power_%s_%04d.pdf'%(car.name, frame)
    if len(glob.glob(outname) ) > 0:
        print "PUNT ON", outname
        continue
    print "GOING TO DO", outname
    all_the_spectra(car,fields)
    sys.exit(0)
sys.exit(1)
