
import equilibriumdata
import numpy as na
import pdb
#equilibriumdata.NumberDensities
#equilibriumdata.GasPressures
akNumberDensity = 1 # cm^-3
pc_in_cm = 3.085677581e18 #cm
Myr_in_s = 3.15576e13 #Myr in s.  Assuming 365.25 days in a year.
akPressure = 7.3014e-13 #g cm^-1 s^-2
akLength = 6.1714e20 #cm
akTime = 9.3e14 #s
mp = 1.67262171e-24  #g, from phys_constants.h
kb = 1.3806504e-16 #erg/K from phys_constants.h 
akTemperature =mp*(akLength/akTime)**2/kb 
time_one_pc_one_kelvin = (kb/mp)**(-1./2)*pc_in_cm #10.78 Myr.
Tunits_pc_myr = mp/kb*(pc_in_cm/Myr_in_s)**2  #115 Kelvin.


DensityUnits = 1.67262171e-24 #1 proton per cc
TimeUnits = 3.15576e14        #Myr in s
LengthUnits = 6.171355162e20  #200 pc in cm
TempKelvin = 1e4
TempUnits = mp*(LengthUnits/TimeUnits)**2/kb 
TempCode =TempKelvin/TempUnits
print "TempCode", TempCode
#print Tunits_pc_myr 
def n_p(n):
    """takes number density *n* (physical units, cm^-3) and returns [p_code,p_phys, p/k phys]"""
    p_phys = na.interp(n,equilibriumdata.NumberDensities,equilibriumdata.GasPressures)
    return [p_phys, p_phys/akPressure, p_phys/kb]
    

def Lam(n,p_phys):
    """Takes number density *n* [cm^-3]and presssure in physical units *p_phys* and returns
    the cooling rate [Lambda, Lamba_code]"""
    Gam = 2e-26 #Gram Centimeter^2/Second^3;
    T = p_phys/(n*kb)
    print "T=",T,p_phys/kb, n
    Lambda= -Gam *(1e7*na.exp((-1.184e5)/(T + 1000)) + 14e-3*na.sqrt(T)*na.exp(-92/T));
    return [Lambda, Lambda*akTime/(akTemperature*kb)]
def dedt(n,p_phys):
    Lambda, Lambda_code = Lam(n,p_phys)
    Gam = 2e-26 #Gram Centimeter^2/Second^3;
    dedt_phys = Lambda*n/mp + Gam/mp
    return [dedt_phys, dedt_phys*akTime/(akTemperature*kb)*mp]

