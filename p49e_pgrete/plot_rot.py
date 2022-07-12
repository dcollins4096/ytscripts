from go import *
import p42_pg
reload(p42_pg)
if './p42_new_driving_ak' not in sys.path:
    sys.path += ['./p42_new_driving_ak']

#car = taxi.load('za01')
car = taxi.load('pg0')
p42_pg.plot_dil(car,100,fieldname='avel')
