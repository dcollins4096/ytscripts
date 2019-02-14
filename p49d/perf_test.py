
from go import *
import performance_tools as pt
reload(pt)

car = taxi.load('p49d_a04')
filename = "%s/performance.out"%car.directory
p = pt.perform(filename)
print(p.data['Stochastic']['Mean Time']/p.data['MHDRK2']['Mean Time'])
