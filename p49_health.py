
#ac19_M1_Ma0.3_512
#ac22_M3_Ma03_512
#ac23_Ma6_Ma01_256
#ac25_M01_Ma01_512
#ac26_M01_Ma5_512

if 'flt' not in dir():
    flt = taxi.fleet(['ac19','ac22','ac23','ac25','ac26'])
    flt('car.qb_load()')
    flt['frames']='last'
    flt('car.get_frames()')
class dumb():
    def __init__(self):
        pass
def print_status(car):
    last_frame=max(car.frame_dict.keys())
    lfd = car.frame_dict[last_frame]
    last_quan = max(car.qb.stuff['frames'])
    if 'EB' in car.qb.stuff:
        last_EB = max(car.qb.stuff['EB'].keys())
    else:
        last_EB = -1
    #{'SetNumber': 242, 'dsname': './DD0242/data0242', 'DirPrefix': './DD', 'time': 1.209999951012777, 'name_files': 'data', 'cycle': 69519}

    print "{c.name}, {lfd[SetNumber]}, {lfd[cycle]}, {lfd[time]}, {:d}, {:d}".format( last_quan, last_EB, c=car, lfd=lfd)
    return lfd
#lfd=print_status(flt[0])
for car in flt:
    print_status(car)
