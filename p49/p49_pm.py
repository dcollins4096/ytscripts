if 'flt' not in dir():
    flt = taxi.fleet(['aa19', 'aa20', 'aa21', 'aa22', 'axb19', 'axb20', 'axb21', 'axb22', 'az19', 'az20', 'az21', 'az22'])
    flt['frames']='last'
    flt('car.load()')
flt['clobber_plot'] = False
flt.save()
#flt('print car.name, max(car.qb.stuff["EB"].keys())')
#flt_x = taxi.fleet(['axb19', 'axb20', 'axb21', 'axb22'])
#flt_a = taxi.fleet(['aa19', 'aa20', 'aa21', 'aa22'])
#flt_z = taxi.fleet(['az19', 'az20', 'az21', 'az22'])
#flt_a['axis'] = [1,2]
#flt_z['axis'] = [1,2]
#flt_a.save()
#flt_z.save()
#for car in flt_x:
#    #car = taxi.taxi(car_name)
#    car.frames='every 10'
#    if len(sys.argv) > 2:
#        car.load()
#        s=nar(sorted(car.frame_dict.keys()))
#        last = s[-1]
#        first = {'ax19':90,'ax20':120,'ax21':110,'ax22':60}[car.name]
#        car.frames = s[slice(first,last,5)].tolist()+[last]
#flt_x['frames']
