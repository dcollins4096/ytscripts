from go import *

if 'always_do' not in dir():
    always_do = False
if 'car' not in dir() or always_do:
    car=taxi.load('p49d_a01')
    car.frames='all'
    frame=0

    ds = car.load(frame)
#car.plot()
quan_list=['time','density_avg','density_std','bx_avg', 'by_avg', 'bz_avg', 'bx_std', 'by_std', 'bz_std', 'vx_avg',  'vy_avg', 'vz_avg',#'vx_std','vy_std',,'density_avg', 'density_std'
          'px_avg','py_avg', 'pz_avg', 'ex_avg', 'ey_avg', 'ez_avg', 'time']#, 'volume',   'vz_std']
qb_list=['t', 'density','density2', 'Bx', 'By', 'Bz','bx2', 'by2', 'bz2','vx', 'vy', 'vz', 'px', 'py', 'pz', 'ex', 'ey', 'ez', 't']# 'Bfield_strength', 'AlfMach', 'beta', 'AlfvenSpeed', 'frames', 'ke_tot', 'ke_rel', 'grav_pot', 'grav_pot_2', 'gas_work', 'tdyn']
quan_list=['time','density_avg'] #,'density_std','bx_avg', 'by_avg', 'bz_avg', 'bx_std', 'by_std', 'bz_std', 'vx_avg',  'vy_avg', 'vz_avg',#'vx_std','vy_std',,'density_avg', 'density_std'
quan_to_qb_map=dict(zip(quan_list,qb_list))

special_list=[  'mach', 'cycle']
from collections import defaultdict
if 'all_quan' not in dir() or always_do:
    print("wtf")
    all_quan = defaultdict(lambda: list())
    dumb=[]
    for n in car.frame_dict:
        this_name = "%s/%s.AverageQuantities.h5"%(car.directory,car.frame_dict[n]['dsname'])
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
import turb_quan
reload(turb_quan)
print("density", all_quan['density_avg'])
if 'qb' not in dir() or always_do:
    qb = turb_quan.quan_box(car)
    qb.plot_format='png'
    qb()
    #qb.plot()

def in_ish(v,arr):
    return (np.abs( v - arr) < 1e-9).any()
if 1:
    slice_qb = np.zeros_like(qb.stuff['t'],dtype='bool')
    slice_quan = np.zeros_like(all_quan['time'],dtype='bool')
    for n, v in enumerate(qb.stuff['t']):
        slice_qb[n] = in_ish(v , all_quan['time']) 
    slice_qb[0]=False
    for n, v in enumerate(all_quan['time']):
        slice_quan[n] = in_ish(v , qb.stuff['t'])
    #print(qb.stuff['t'][slice_qb])
    #print(all_quan['time'][slice_quan])
    #print(qb.stuff['t'][slice_qb]- all_quan['time'][slice_quan])

    for q in quan_list:
        name_qb =  quan_to_qb_map[q]
        qbv =  nar(qb.stuff[ name_qb])[slice_qb] 
        qv =  all_quan[q][slice_quan]
        if name_qb in ['density2']: #, 'bx2','by2','bz2']:
            qbv = np.sqrt(qbv)
        print("=== Q %s (%s) sum %0.2e diff %0.2e==="%(q, name_qb, qv.sum(), np.abs(qbv-qv).sum()))
        if name_qb in ['bx2']:
            print(qbv)
            print(qv)
            print( qbv -qv)

