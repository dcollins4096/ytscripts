from go import *

if 'car' not in dir():
    car=taxi.taxi('p49d_a02')
    car.frames='all'
    frame=3
    ds = car.load(frame)
#car.plot()
quan_list=['time','density_avg','density_std','bx_avg', 'by_avg', 'bz_avg', 'bx_std', 'by_std', 'bz_std', 'vx_avg',  'vy_avg', 'vz_avg',#'vx_std','vy_std',,'density_avg', 'density_std'
          'px_avg','py_avg', 'pz_avg', 'ex_avg', 'ey_avg', 'ez_avg', 'time']#, 'volume',   'vz_std']
qb_list=['t', 'density','density2', 'Bx', 'By', 'Bz','bx2', 'by2', 'bz2','vx', 'vy', 'vz', 'px', 'py', 'pz', 'ex', 'ey', 'ez', 't']# 'Bfield_strength', 'AlfMach', 'beta', 'AlfvenSpeed', 'frames', 'ke_tot', 'ke_rel', 'grav_pot', 'grav_pot_2', 'gas_work', 'tdyn']
quan_to_qb_map=dict(zip(quan_list,qb_list))

special_list=[  'mach', 'cycle']
from collections import defaultdict
if 'all_quan' not in dir() or True:
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
import turb_quan
reload(turb_quan)
if 'qb' not in dir():
    qb = turb_quan.quan_box(car)
    qb.plot_format='png'
    qb()
    qb.plot()

if 1:
    slice_quan = slice(None)
    slice_qb = slice(1,None)
    print(qb.stuff['t'][slice_qb])
    print(all_quan['time'][slice_quan])
    print(qb.stuff['t'][slice_qb]- all_quan['time'][slice_quan])
    for q in quan_list:
        name_qb =  quan_to_qb_map[q]
        qbv =  qb.stuff[ name_qb][slice_qb] 
        qv =  all_quan[q][slice_quan]
        if name_qb in ['density2']: #, 'bx2','by2','bz2']:
            qbv = np.sqrt(qbv)
        print("=== Q %s (%s) sum %0.2e diff %0.2e==="%(q, name_qb, qv.sum(), np.abs(qbv-qv).sum()/np.mean(qv)))
        if name_qb in ['bx2']:
            print(qbv)
            print(qv)
            print( qbv -qv)

