

data = at.leaf_clumps[3][0].data
AvgVelocity = data.ds.arr([(data['cell_mass']*data[v]).sum()/data['cell_mass'].sum()
                              for v in ['x-velocity','y-velocity','z-velocity']])
data.set_field_parameter("bulk_velocity",AvgVelocity)
print data['VelocityDispersionSquared']
