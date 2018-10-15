
def scaling(ds,n_cgs,T_K=10):
    """n particles/cm^3
    T_K in Kelvin, defaults to 0.2km/s.  Actually totally unused"""
    FourPi_G_code = ds['GravitationalConstant']
    FourPi_G_CGS = 4*np.pi*6.7e-8 # cm^3/g/s^2
    sound_speed_cgs = 0.2 #km/s
    mp = 1.67e-24 #gram
    mu = 2.3     
    rho_cgs = mu*mp*n_cgs
    rho_code = 1.0
    MachNumber = 9 #i think
    to_seconds = (FourPi_G_code*rho_code/(FourPi_G_CGS*rho_cgs))**0.5
    pc_to_km = 3.0e13 #km/pc
    v_rms_cgs = MachNumber*sound_speed_cgs
    code_length = v_rms_cgs * to_seconds/pc_to_km
    print code_length

scaling(ds_late,1000)
