execfile('go')
ef('turb_quan.py')
import taxi
car = taxi.taxi(sys.argv[1])
frames=None
if len(sys.argv) > 2:
    if sys.argv[2] == 'car':
        pass #just leave the frames as they are on disk.
    else:
        try:
            #I don't really care for using try-except for this, but 
            #I can't think of a better way...
            car.frames = [int(a) for a in sys.argv[2:]]
        except:
            L = len(sys.argv[2:])
            car.frames="%s "*L%tuple(sys.argv[2:])
            car.frames=car.frames[:-1]
print "FRAMES", car.frames
print "FRAMES", car.return_frames()
#car.plot(0)
if 1:
    ef('turb_quan.py')
    qb = quan_box(car)
    qb.plot_format='png'
    qb.load()
    qb.EBall()
    #qb()
    #qb.plot()
    #pickle_name = 'quan_box_%s.pickle'%car.name
    #extant_quan=None
    #if len(glob.glob(pickle_name)):
    #    print "LOAD", pickle_name
    #    extant_quan = fPickle.load(pickle_name)
    #car.axis=['x','y','z']
    #car.fields=['density']
    #car.plot()
