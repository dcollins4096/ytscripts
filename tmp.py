
if 'ef' not in dir():
    execfile('go')
reload(taxi)
import turb_quan
reload(turb_quan)

import xtra_energy_fields
reload(xtra_energy_fields)
car = taxi.taxi('ab25')
car.load()
field = 'eng_x'
#car.ds.add_field(field,**xtra_energy_fields.field_args[field])
xtra_energy_fields.dave_add_field(car.ds) #, field_name = field)
print(car.ds.index.grids[0][field])
#qb = turb_quan.quan_box(car)
#qb.load()
field_args={}
