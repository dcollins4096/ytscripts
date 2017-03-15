from xtra_operators import *    
#yt.add_field('scaled_div_b',  function=_scaled_div_b, validators=[yt.ValidateGridType()])

# \dU/dt = -grad P - 1/8 pi \grad(b^2) -1/4pi B\cdot \grad 
def _od(field,data):
    return data['density'].in_units('code_density').v
yt.add_field('od',_od)









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
