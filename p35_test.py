
if 'flt' not in dir():
    flt = taxi.fleet(['n01','n03'])
fields = ["Dark_Matter_Density", "Density", "Electron_Density", "GasEnergy", "H2II_Density",
"H2I_Density", "HII_Density", "HI_Density", "HM_Density", "HeIII_Density", "HeII_Density", "HeI_Density", "Temperature", "TotalEnergy"]

flt['frames'] = range(27)
flt['Colorbar']='fixed'
for field in fields[0:]:
    #flt.find_extrema([field])
    #flt.
    flt.phase(['density',field,'cell_mass'])
