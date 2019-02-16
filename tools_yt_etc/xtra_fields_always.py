
#yt.add_field('scaled_div_b',  function=_scaled_div_b, validators=[yt.ValidateGridType()])


# \dU/dt = -grad P - 1/8 pi \grad(b^2) -1/4pi B\cdot \grad 
def _od(field,data):
    return data['density'].in_units('code_density').v
yt.add_field('od',_od)

if 0:
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
