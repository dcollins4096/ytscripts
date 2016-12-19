execfile('go')
execfile('p42_pg.py')
execfile('p42_spectra.py')
ef('p42_quan.py')
try:
    import taxi
    car = taxi.taxi(sys.argv[1])
#car.plot(0)
    if 1:
        car.fill(0)
        all_frames = car.frame_dict.keys()
        car.frames = all_frames[::10]
        if all_frames[-1] not in car.frames:
            car.frames += [all_frames[-1]]
        car.axis=['x','y','z']
        car.fields=['density']
        car.plot()
        qb = quan_box(car)
        qb(car)
        fPickle(qb,'quan_box_%s.pickle'%car.name)
except:
    raise
finally:
    fptr = open("p49_batch_%s_finished"%car.name)
    fptr.write("finished\n")
    fptr.close()
