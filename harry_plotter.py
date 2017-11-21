execfile('go')
reload(taxi)
if 'car_name' not in dir():
    car_name = sys.argv[1]
car = taxi.taxi(car_name)
car.clobber_plot=False
car.frames='every 10'
if len(sys.argv) > 2:
    car.load()
    s=nar(sorted(car.frame_dict.keys()))
    last = s[-1]
    first = {'ax19':90,'ax20':120,'ax21':300,'ax22':60}[car.name]
    car.frames = s[slice(first,last,5)].tolist()+[last]
print "XXXXXXX", car.return_frames()
car.plot()

