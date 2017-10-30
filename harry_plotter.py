execfile('go')
reload(taxi)
if 'car_name' not in dir():
    car_name = sys.argv[1]
car = taxi.taxi(car_name)
if len(sys.argv) > 2:
    car = taxi.taxi(car_name)
    car.load()
    s=nar(sorted(car.frame_dict.keys()))
    last = s[-1]
    first = {'ax19':80,'ax20':120,'ax21':110,'ax22':55}[car.name]
    car.frames = s[slice(first,last,5)].tolist()+[last]
print car.return_frames()
car.plot()

