from xtra_operators import *    
eng_units = 'g/(cm*s**2)'

ef('xtra_fields_always.py')
if 'add_field' not in dir():
    from yt import add_particle_filter, add_field, ValidateGridType

def formed_star(pfilter, data):
    filter = data["all", "creation_time"] > 0
    return filter

add_particle_filter("formed_star", function=formed_star, filtered_type='all',
                    requires=["creation_time"])

ef('xtra_fields_particles.py')






if 0:
    #these take arguments (CenterOfMass, AngularMomentumVector, Positions, Velocities)
    import  yt.utilities.math_utils as mu
    compute_vphi = mu.compute_rotational_velocity 
    compute_vz = mu.compute_parallel_velocity 
    compute_vr = mu.compute_radial_velocity 
    #'compute_parallel_velocity', 'compute_radial_velocity', 'compute_rotational_velocity'

    def set_cylindrical_stuff(data, COM=None, L = None):
        if COM is not None:
            CoM = data.quantities['CenterOfMass']()
        if L is not None:
            L   = data.quantities['AngularMomentumVector']()
        data.set_field_parameter('CoM',CoM)
        data.set_field_parameter('L',L)
        
    def _vphi_math_util(field,data):
        """doesn't work well."""
        if not data.has_field_parameter('CoM'):
            out = np.zeros_like(data['x-velocity'])
            return out
        CoM = data.get_field_parameter('CoM')
        if type(CoM) is types.FloatType:
            out=np.zeros_like(data['x-velocity'])
            return out
        L   = data.get_field_parameter('L')
        P = nar(zip( data['x'], data['y'], data['z'] ))
        V = nar(zip( data['velocity_x'], data['velocity_y'], data['velocity_z']))
        out = compute_vphi(CoM,L, P,V)
        return out
        
    def _vphi(field,data):
        """doesn't work well."""
        if not data.has_field_parameter('CoM'):
            out = np.zeros_like(data['x-velocity'])
            return out
        CoM = data.get_field_parameter('CoM')
        if type(CoM) is types.FloatType:
            out=np.zeros_like(data['x-velocity'])
            return out
        L   = data.get_field_parameter('L')
        P = nar(zip( data['x'], data['y'], data['z'] ))
        V = nar(zip( data['velocity_x'], data['velocity_y'], data['velocity_z']))
        out = compute_vphi(CoM,L, P,V)
        return out
    yt.add_field('vphi',function=_vphi, validators=[yt.ValidateParameter(v) for v in ['CoM','L']])


#print "MONKEY ON THE CAR"
