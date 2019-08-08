from go import *
import p49_fields
import turb_quan
def EBshorty(car):
    if type(car) is str:
        car = taxi.load(car)
    car.derived_fields['QU'] = p49_fields.add_QU
    car.fields=['density','magnetic_field_strength']
    car.axis=[0]
    car.frames='last'
    #car.plot()
    qb = turb_quan.quan_box(car=car)
    qb.EBall()

flt = taxi.fleet([ "za01", "zb01", "zc01_quan", "zd01_quan", "ze01_quan"][::-1])
for car in flt:
    EBshorty(car)

