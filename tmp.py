if 'ef' not in dir():
    execfile('go')
    
density = 3.85e-21
vel = 1.8e4
field_unit="%0.2e Gauss"%(np.sqrt(4*np.pi*density*vel**2))
vel_unit = '%0.2e cm/s'%vel
units={'Density':"%0.2e g/cm^3"%3.85e-21,'x-velocity':vel_unit,'y-velocity':vel_unit,'z-velocity':vel_unit,
       'Bx':field_unit,'By':field_unit,'Bz':field_unit}
