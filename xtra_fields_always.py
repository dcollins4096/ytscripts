
#yt.add_field('scaled_div_b',  function=_scaled_div_b, validators=[yt.ValidateGridType()])


# \dU/dt = -grad P - 1/8 pi \grad(b^2) -1/4pi B\cdot \grad 
def _od(field,data):
    return data['density'].in_units('code_density').v
yt.add_field('od',_od)


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
yt.add_field("timestep",function=dt_enzo, validators=[yt.ValidateGridType()]) #, units='dimensionless')

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
