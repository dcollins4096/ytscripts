
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

def doit(flt,fun):
    for car in flt:
        fun(car)

def check_output(car):
    last_frame=max(car.frame_dict.keys())
    last_set = car.frame_dict[last_frame]['dsname']
    last_set_full = "%s/%s"%(car.directory,last_set)
    print last_set_full
    print car.name, glob.glob(last_set_full)

doit(flt, check_output)

def check_quan(car, print_ok=True, print_not_ok = False):
   file_list = glob.glob("%s/DD????"%car.directory)
   for f in file_list:
       frame = int(f[-4:])
       if frame in car.qb.stuff['frames'] and frame in car.qb.stuff['EB']:
           if print_ok:
               a,b=f.split('/')[-2:]
               print "%s/%s.tar"%(a,b)
       else:
           if print_not_ok: 
               print "NOT OK"



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
def update(car):
    new_pickle_name = 'xfer/quan_box_%s.pickle'%car.name
    if glob.glob(new_pickle_name) != []:
        car.qb.stuff.update(fPickle.load(new_pickle_name))
    else:
        print "no new pickle found"
if 'do_update' in dir():
    for car in flt:
        print_status(car)
        update(car)
        print_status(car)

    if 'yes_save' in dir():
        for car in flt:
            car.qb.dump()
            print "dumped", car.name

