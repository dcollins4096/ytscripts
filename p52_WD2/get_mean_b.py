
from go import *
reload(taxi)
cars = ['p52_441','p52_432','p52_433','p52_434']
if 'fltb' not in dir():
    fltb = taxi.fleet(cars)
labelmap = {'p52_441':'D22','p52_432':'D26','p52_433':'D27','p52_434':'D28'}


means={}
std={}
for car in fltb.taxi_list:
    ds = car.load(0)
    cg=ds.all_data()
    Bx=cg['Bx'].v
    By=cg['By'].v
    Bz=cg['Bz'].v
    Bx0 = Bx.mean()
    By0 = By.mean()
    Bz0 = Bz.mean()
    means[car.name]=nar([Bx0,By0,Bz0])
    std[car.name] = nar([Bx.std(), By.std(), Bz.std()])

for car in fltb.taxi_list:
    Bx0,By0,Bz0=means[car.name]

    print( "===")
    print( "%0.2e %0.2e %0.2e"%(Bx0,By0,Bz0))
    print( "%0.2e %0.2e %0.2e"%tuple( means[car.name]/std[car.name]))

