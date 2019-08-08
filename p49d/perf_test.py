
from go import *
import performance_tools as pt
reload(pt)

if 0:
    car = taxi.load('p49d_a04')
    filename = "%s/performance.out"%car.directory
    p = pt.perform(filename)
#print(p.data['Stochastic']['Mean Time']/p.data['MHDRK2']['Mean Time'])

if 'p1' not in dir():
    car_base = 'ze01'
    c1 = taxi.load(car_base)
    c2 = taxi.load(car_base+"_quan")
    filename1 = "%s/performance.out"%c1.directory
    filename2 = "%s/performance.out"%c2.directory
    p1 = pt.perform(filename1)
    p2 = pt.perform(filename2)


smooth_len = 11
plt.clf()
y1 = pt.smooth( p1.data['Total']['Mean Time'],smooth_len)
y2 = pt.smooth( p2.data['Total']['Mean Time'],smooth_len)
y3 = pt.smooth( p2.data['AverageQuantities']['Mean Time'],smooth_len)
plt.plot(p1.data['Total']['Cycle'],y1,label=car_base)
plt.plot(p2.data['Total']['Cycle'],y2,label=car_base+"_quan")
plt.plot(p2.data['Total']['Cycle'],y1[:len(y3)]+y3,label=car_base+"_quan est")
plt.legend(loc=0)
plt.ylim(0,11)
outname = "P49d_perf_total_%s.png"%car_base
plt.savefig(outname)
print("plot "+outname)
