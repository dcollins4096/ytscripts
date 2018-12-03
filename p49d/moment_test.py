
if 'car' not in dir():
    car=taxi.taxi('p49d_a02')
#car.plot()
if 'frame' not in dir():
    frame=3
ds = car.load(frame)
ad = ds.all_data()
dv = ad['cell_volume']
vtotal = dv.sum()
field='Bx'
mean_field=(ad[field]*dv).sum()/vtotal
var_field=( (ad[field]-mean_field)**2*dv).sum()/vtotal
sqr_field=( (ad[field])**2*dv).sum()/vtotal
std_field = np.sqrt(var_field)
print("Mean field %0.7e"%(mean_field))
print("Msq  field %0.7e"%(mean_field**2))
print("var  field %0.7e"%(var_field))
print("std  field %0.7e"%(std_field))
print("sqr  field %0.7e"%(sqr_field))
print("std  field %0.7e"%( np.sqrt(sqr_field-mean_field**2)))
reg = ad
total_volume = reg['cell_volume'].sum()
#print((( reg[field].in_units('code_field').v)**2*reg['cell_volume']).sum()/total_volume)  
#print((( reg[field].in_units('code_field').v-mean_field.in_units('code_field').v)**2*reg['cell_volume']).sum()/total_volume)  
print(np.sqrt((( ad[field]-mean_field)**2*dv).sum()/vtotal)  )
if 0:
    mean_density=(ad['density']*dv).sum()/vtotal
    var_density=( (ad['density']-mean_density)**2*dv).sum()/vtotal
    sqr_density=( (ad['density'])**2*dv).sum()/vtotal
    std_density = np.sqrt(var_density)
    print("Mean density %0.7e"%(mean_density))
    print("Msq  density %0.7e"%(mean_density**2))
    print("var  density %0.7e"%(var_density))
    print("std  density %0.7e"%(std_density))
    print("sqr  density %0.7e"%(sqr_density))
    print("std  density %0.7e"%( np.sqrt(sqr_density-mean_density**2)))
    reg = ad
    total_volume = reg['cell_volume'].sum()
    print((( reg['density'].in_units('code_density').v)**2*reg['cell_volume']).sum()/total_volume)  
    print((( reg['density'].in_units('code_density').v-mean_density.in_units('code_density').v)**2*reg['cell_volume']).sum()/total_volume)  
    print((( ad['density']-mean_density)**2*dv).sum()/vtotal)  
