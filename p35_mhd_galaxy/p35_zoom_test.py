reload(taxi)
fleet=taxi.fleet(['n05','n06','n08'])
fleet['frames']=[12]
fleet('car.set_center_max()')
fleet['fields']=['Temperature']
fleet['weight_field']='density'
fleet['zoom_sequence']=[(n,'kpc') for n in [15,10,5,1,0.5,0.1,0.05]]
fleet.plot()


