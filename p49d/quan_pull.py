
if 'car' not in dir():
    car=taxi.taxi('p49d_a02')
#car.plot()
if 'frame' not in dir():
    frame=3
ds = car.load(frame)

from collections import defaultdict
all_quan = defaultdict(lambda: list())
for n in car.frame_dict:
    this_name = "%s/%s.AverageQuantities.h5"%(car.directory,car.frame_dict[n]['dsname'])
    if not os.path.exists(this_name):
        print("No such file, skipping: %s"%this_name)
        continue
    fptr = h5py.File(this_name,'r')
    for key in fptr:
        all_quan[key].append( fptr[key][:])
    try:
        pass
    except:
        raise
    finally:
        fptr.close()
for key in all_quan:
    all_quan[key] = np.array( all_quan[key]).flatten()

