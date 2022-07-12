from go import *
print("WTF MATE")
def dt_picker( car, target_times):
    fptr = open( "%s/OutputLog"%car.directory,'r')
    lines=fptr.readlines()
    fptr.close
    names = []
    times = []
    for line in lines:
        sp = line.split()
        names.append(sp[2])
        times.append(sp[4])
    times = np.array(times,dtype='float')
    out_frames = []
    out_names = []
    out_times = []
    for t in target_times:
        if t < times.min():
            continue
        if t > times.max():
            continue
        dist = np.abs( t - times) 
        a = np.where( dist ==dist.min())[0][0]
        try:
            out_frames.append(a)
            out_names.append( names[a])
            out_times.append(times[a])
        except:
            pdb.set_trace()

    return out_frames,out_names,out_times

car = taxi.load('ze01')
out=dt_picker(car,target_times=np.arange(0,0.4,0.01))
