from go import *
import xtra_dynamo_fields as xdf
reload(xdf)
import xtra_energy_fields as xe
reload(xe)


cars = ['p52_441','p52_432','p52_433','p52_434']
if 'flt' not in dir():
    flt = taxi.fleet(cars)
labelmap = {'p52_441':'D22','p52_432':'D26','p52_433':'D27','p52_434':'D28'}

flt['derived_fields'] = {'energy_terms':xe.add_extra_energies,'dynamo_terms':xdf.add_dynamo_fields }
flt['fields'] = ['Ltension_squared']
flt['name_syntax'] = 'preset'
flt['name_dir'] = 'datasets/EE'
flt['name_files'] = 'data'
flt['frames'] = [0]
if 0:
    for car in flt.taxi_list[0:1]:
        ds=car.load(0)
        ad=ds.all_data()
        print(np.abs(ad['Ltension_squared']).sum())
flt['operation'] = 'CenterSlice'
flt['callbacks'] = ['time_title']
flt['frames'] = [0,10,60]
flt.plot()
    

