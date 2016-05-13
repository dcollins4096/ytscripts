def _Rel_KineticEnergy(field,data):
    bv_x,bv_y,bv_z = data.get_field_parameter("bulk_velocity")
    return 0.5 * (data["Density"] * (
        (data["x-velocity"] - bv_x)**2
        + (data["y-velocity"] - bv_y)**2
        + (data["z-velocity"] - bv_z)**2 )).v

yt.add_field("RelKineticEnergy",function = _Rel_KineticEnergy, 
                validators=[yt.ValidateParameter("bulk_velocity")])

def _MassToFluxNonCritical(field,data):
    MassToFlux = data['cell_mass']/(data['magnetic_field_strength']*(data['dx'])**2)
    return MassToFlux.v
yt.add_field('MassToFluxNonCritical',function = _MassToFluxNonCritical)



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

yt.add_field("VelocityDispersionSquared", function=_VelocityDispersionSquared,units='code_velocity**2',
                validators=[yt.ValidateParameter("bulk_velocity")], take_log = False)



#print "MONKEY ON THE CAR"
