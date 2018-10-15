
from yt.utilities import physical_constants as units
#if 'fa01' not in dir()
fa01=taxi.taxi('fa01')
fa02=taxi.taxi('fa02')

car = fa02
car.name_syntax='preset'
car.name_files='data'
car.frames=range(0,150,10)
def time_series(car,field):
    print car.frames
    min_series=[]
    max_series=[]
    avg_series=[]
    for frame in car.frames:
        reg=car.get_region(frame)
        if 0:
            fmin, fmax = reg.quantities['Extrema'](field)
            if np.abs(1-fmin/fmax) > 1e-6:
                favg = reg.quantities['WeightedAverageQuantity'](field,'cell_volume')
            else:
                favg = 0.5*(fmin+fmax)
        else:
            fmin = reg[field][ reg[field].size/2]
            fmax=fmin
            favg=fmin
        min_series.append(fmin)
        max_series.append(fmax)
        avg_series.append(favg)
    return nar(min_series), nar(max_series), nar(avg_series)

if 0:
    min_pressure, max_pressure, avg_pressure = time_series(car,'pressure')
    min_T, max_T, avg_T = time_series(car,'Temperature')
    min_density, max_density, avg_density = time_series(car,'density')
    len(min_density)

print min_T
import equilibriumdata
NumberDensities = nar(equilibriumdata.NumberDensities)
GasPressure = nar(equilibriumdata.GasPressures)

MeanAtomicMass = 0.7
fields=['density','GasPressure']
if 'davesNumberDensity' in fields:
	Densities = NumberDensities
elif 'density' in fields:
	Densities = NumberDensities *MeanAtomicMass*units.mp

if 'GasPressure' in fields:
	Pressures = GasPressure
elif 'PoverK' in fields:
	Pressures = GasPressure/units.kboltz

plt.clf()
def PfromTempDen(Temp,Den):
    mp = 1.67262171e-24  #g, from phys_constants.h
    kb = 1.3806504e-16 #erg/K from phys_constants.h 
    return Den*Temp*kb/mp
def TfromPempDen(P,Den):
    mp = 1.67262171e-24  #g, from phys_constants.h
    kb = 1.3806504e-16 #erg/K from phys_constants.h 
    return P/(Den*kb/mp)
plt.plot(Densities,Pressures,c='b')
plt.plot(Densities,PfromTempDen(1e4,Densities),c='k')
plt.plot(Densities,PfromTempDen(5e3,Densities),c='r')
plt.plot(Densities,PfromTempDen(1e3,Densities),c='k')
plt.plot(Densities,PfromTempDen(1e2,Densities),c='k')
rmap = rainbow_map(len(min_density))
c=[rmap(i) for i in range(len(min_density))]
sct=plt.scatter(min_density,min_pressure, c=c)
plt.xlabel('density')
plt.ylabel('pressure [erg/cm^3]')
plt.xscale('log'); plt.yscale('log')
plt.ylim(5e-14)
outname = "p6_test_%s.pdf"%car.outname
plt.savefig(outname)
print outname 
