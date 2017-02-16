execfile('go')
ef('turb_quan.py')
import taxi
car = taxi.taxi(sys.argv[1])
frames=None
if len(sys.argv) > 2:
    if sys.argv[2] == 'car':
        frames = car.frames
    else:
        frames = [int(a) for a in sys.argv[2:]]
    
#car.plot(0)
if 1:
    pickle_name = 'quan_box_%s.pickle'%car.name
    extant_quan=None
    if len(glob.glob(pickle_name)):
        print "LOAD", pickle_name
        extant_quan = fPickle.load(pickle_name)
    if frames is None:
        car.fill(0)
        all_frames = car.frame_dict.keys()
        car.frames = all_frames[::10]
        if all_frames[-1] not in car.frames:
            car.frames += [all_frames[-1]]
    else:
        car.frames=frames

    #car.axis=['x','y','z']
    #car.fields=['density']
    #car.plot()
    qb = quan_box(car)
    qb(car,extant_quan=extant_quan)
    print "SAVE", pickle_name
    fPickle.dump(qb,pickle_name)
