from go import *
car = taxi.load('p52_b06')

car.operation='CenterSlice'
car.frames='all'
car.fields=['density','pressure','velocity_magnitude']
#car.outname = "P52_ic/%s"%car.name
car.plot()
