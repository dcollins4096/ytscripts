
#yt.add_field('scaled_div_b',  function=_scaled_div_b, validators=[yt.ValidateGridType()])


# \dU/dt = -grad P - 1/8 pi \grad(b^2) -1/4pi B\cdot \grad 
def _od(field,data):
    return data['density'].in_units('code_density').v
yt.add_field('od',_od)

def dt_rk2(field,data):
    """Casting all components to code units eliminates some source of reduced precision that
    was being incurred during thing 1/dx step (I think?)"""
    aye=1.0
    if data.ds['ComovingCoordinates']:
        aye=(data.ds['CosmologyInitialRedshift']+1)/\
            ( data.ds['CosmologyCurrentRedshift']+1)

    if hasattr(data.dds,'in_units'):
        dx, dy, dz = data.dds.in_units('code_length').v
    else:
        dx, dy, dz = data.dds
    dxinv = 1./(aye*dx)
    dyinv = 1./(aye*dy)
    dzinv = 1./(aye*dz)
    cs2 = (data['sound_speed'].in_units('code_velocity').v)**2
    Bx = data['Bx']
    By = data['By']
    Bz = data['Bz']
    if hasattr(Bx,'in_units'):
        Bx = Bx.in_units('code_magnetic').v
        By = By.in_units('code_magnetic').v
        Bz = Bz.in_units('code_magnetic').v
    vx = data['velocity_x'].in_units('code_velocity').v
    vy = data['velocity_y'].in_units('code_velocity').v
    vz = data['velocity_z'].in_units('code_velocity').v
    rho = data['density'].in_units('code_density').v
    B2 = Bx*Bx+By*By+Bz*Bz

    temp1 = cs2 + B2/rho
    ca2 = Bx*Bx/rho;
    cf2 = 0.5 * (temp1 + np.sqrt(temp1*temp1 - 4.0*cs2*ca2));
    cf = np.sqrt(cf2);
    v_signal_x = (cf + np.abs(vx));
    GridRank = 3 #fix later
    if GridRank > 1:
      ca2 = By*By/rho;
      cf2 = 0.5 * (temp1 + np.sqrt(temp1*temp1 - 4.0*cs2*ca2));
      cf = np.sqrt(cf2);
      v_signal_y = (cf + np.abs(vy));
    
    if GridRank > 2:
      ca2 = Bz*Bz/rho;
      cf2 = 0.5 * (temp1 + np.sqrt(temp1*temp1 - 4.0*cs2*ca2));
      cf = np.sqrt(cf2);
      v_signal_z = (cf + np.abs(vz));
    
    dt_x = v_signal_x * dxinv;
    dt_y = v_signal_y * dyinv
    dt_z = v_signal_z * dzinv
    
    dt = np.max([dt_x, dt_y, dt_z],axis=0);
    return data.ds.parameters['CourantSafetyNumber'] /dt
#yt.add_field("timestep_rk2",function=dt_rk2, validators=[yt.ValidateGridType()]) #, units='dimensionless')

def dt_enzo(field,data):
    """Casting all components to code units eliminates some source of reduced precision that
    was being incurred during thing 1/dx step (I think?)"""
    Cs = data['sound_speed'].in_units('code_velocity').v
    aye=1.0
    if data.ds['ComovingCoordinates']:
        aye=(data.ds['CosmologyInitialRedshift']+1)/\
            ( data.ds['CosmologyCurrentRedshift']+1)

    if hasattr(data.dds,'in_units'):
        dx, dy, dz = data.dds.in_units('code_length').v
    else:
        dx, dy, dz = data.dds
    #using harmonic mean.
    dti  = (Cs + na.abs(  data['velocity_x'].in_units('code_velocity').v )) /dx
    dti  += (Cs + na.abs( data['velocity_y'].in_units('code_velocity').v ))/dy
    dti  += (Cs + na.abs( data['velocity_z'].in_units('code_velocity').v ))/dz
    dti /= data.ds.parameters['CourantSafetyNumber'] #yes divided: still in reciporical space
    return aye/dti
#yt.add_field("timestep",function=dt_enzo, validators=[yt.ValidateGridType()]) #, units='dimensionless')

def _VelocityDispersionSquared(field, data):
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
        print "No Bulk Velocity"
    #print "dispersion bulk velocity, eh:", bulk_velocity
    new_field =( (data["x-velocity"]-bulk_velocity[0])**2 + 
                 (data["y-velocity"]-bulk_velocity[1])**2 + 
                 (data["z-velocity"]-bulk_velocity[2])**2 )
    return new_field

yt.add_field("VelocityDispersionSquared", function=_VelocityDispersionSquared,units=r"cm**2/s**2",
                validators=[yt.ValidateParameter("bulk_velocity")], take_log = False)


def _rel_kinetic_energy(field,data):
    if data.has_field_parameter('bulk_velocity'):
        vbar = data.get_field_parameter('bulk_velocity')
    else:
        vbar = data.ds.arr([0]*3,'code_velocity')


    vx = data['velocity_x']-vbar[0]
    vy = data['velocity_y']-vbar[1]
    vz = data['velocity_z']-vbar[2]
    return 0.5*data['density']*(vx*vx+vy*vy+vz*vz)
yt.add_field('rel_kinetic_energy',function=_rel_kinetic_energy,units=eng_units,validators=[yt.ValidateParameter('bulk_velocity')])
