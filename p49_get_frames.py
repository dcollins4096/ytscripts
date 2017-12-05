

def read_frames(car):
    frame_file = open(car.directory+'/frames','r')
    frames = []
    for line in frame_file:
        frames.append( int(line[:-1]))
    car.frames=frames


if 'car_list' not in dir():
    car_list=[]
for name in ['ab26','ac26']:
    car=taxi.taxi(name)
    read_frames(car)
    print car.frames
    car.save(name)
    car_list.append(car)
