from GL import *
from collections import defaultdict

def return_average_quantities(file_list=[], all_quan = defaultdict(lambda: list())):
    for this_name in file_list:
        if not os.path.exists(this_name):
            print("No such file, skipping: %s"%this_name)
            continue
        fptr = h5py.File(this_name,'r')
        for key in fptr:
            mything =  list(fptr[key][:].flatten())
            all_quan[key] += mything
        try:
            pass
        except:
            raise
        finally:
            fptr.close()
    for key in all_quan:
        all_quan[key] = np.array( all_quan[key]).flatten()
    return all_quan

def all_quan_from_taxi(car):
    dumb=[]
    file_list=[]
    for n in car.frame_dict:
        file_list.append("%s/%s.AverageQuantities.h5"%(car.directory,car.frame_dict[n]['dsname']))
    all_quan = return_average_quantities(file_list)
    return all_quan
