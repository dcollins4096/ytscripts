
"""
all_x=['Bstretching',
'Bcompression',
'Badvection',
'Ltension',
'Lpressure']



print "MONKEY ON THE CAR"

all_units=["norm_unit_velocity_x", "norm_unit_momentum_x", "norm_unit_magnetic_field_x", \
           "norm_unit_velocity_y", "norm_unit_momentum_y", "norm_unit_magnetic_field_y",\
           "norm_unit_velocity_z", "norm_unit_momentum_z", "norm_unit_magnetic_field_z"]
all_mag=["Bstretching_f", "Bcompression_f", "Badvection_f", "Ltension_f", "Lpressure_f",\
         "Bstretching_p", "Bcompression_p", "Badvection_p", "Ltension_p", "Lpressure_p",\
         "Bstretching_m", "Bcompression_m", "Badvection_m", "Ltension_m", "Lpressure_m",\
         "Bstretching_t", "Bcompression_t", "Badvection_t", "Ltension_t", "Lpressure_t"]
"""

    
from xtra_operators import *
import numpy as np
import yt

def add_dynamo_fields(obj):
    add_field = obj.add_field
    ValidateSpatial = yt.ValidateSpatial
    def _InductionNorm(field,data):
        """ || Curl( V cross B) ||"""
        B = [ data['Bx'],
              data['By'],
              data['Bz'] ]
        V = [ data['x-velocity'],
              data['y-velocity'],
              data['z-velocity'] ]
        E = Cross(V,B)
        dx = data['dx'].flat[0]
        dBdT= Curl( E, [dx]*3)
        return np.sqrt( dBdT[0]**2+dBdT[1]**2 +dBdT[2]**2)
    add_field('InductionNorm',function=_InductionNorm,
                    validators=[ValidateSpatial(1,['Bx','By','Bz', 'x-velocity','y-velocity','z-velocity'])])

    std_validators = [ValidateSpatial(1,['Bx','By','Bz', 'x-velocity','y-velocity','z-velocity'])]
    std_validators_2 = [ValidateSpatial(1,['magnetic_field_x','magnetic_field_y','magnetic_field_z', 'x-velocity','y-velocity','z-velocity'])]
    momentum_rate_units = 'g/(cm**2*s**2)' #momentum density per time = (gram cm/s)/cm^3/time = g /(cm^2
    def _Ltension_x(field,data):
        return  AdotDel(data,['magnetic_field_x','magnetic_field_y','magnetic_field_z'],'magnetic_field_x')/(4*np.pi)
    add_field('Ltension_x', function=_Ltension_x,validators=std_validators_2, take_log=False, units = momentum_rate_units,sampling_type='cell')
    def _Ltension_y(field,data):
        return  AdotDel(data,['magnetic_field_x','magnetic_field_y','magnetic_field_z'],'magnetic_field_y')/(4*np.pi)
    add_field('Ltension_y', function=_Ltension_y,validators=std_validators_2, take_log=False, units = momentum_rate_units,sampling_type='cell')
    def _Ltension_z(field,data):
        return  AdotDel(data,['magnetic_field_x','magnetic_field_y','magnetic_field_z'],'magnetic_field_z')/(4*np.pi)
    add_field('Ltension_z', function=_Ltension_z,validators=std_validators_2, take_log=False, units = momentum_rate_units,sampling_type='cell')
    def _Lpressure_x(field,data):
        return -1*grad(data,'magnetic_energy',0) #the 1/8pi is in the magnetic energy.
    add_field('Lpressure_x', function=_Lpressure_x,validators=std_validators, take_log=False, units=momentum_rate_units,sampling_type='cell')
    def _Lpressure_y(field,data):
        return -1*grad(data,'magnetic_energy',1)
    add_field('Lpressure_y', function=_Lpressure_y,validators=std_validators, take_log=False, units=momentum_rate_units,sampling_type='cell')
    def _Lpressure_z(field,data):
        return -1*grad(data,'magnetic_energy',2)
    add_field('Lpressure_z', function=_Lpressure_z,validators=std_validators, take_log=False, units=momentum_rate_units,sampling_type='cell')

    magnetic_rate_unit = 'code_magnetic/s'
#On these units:  Using Bx and magnetic_field_x are interchangeable; the units ('code_magnetic' vs 'gauss') will take care of 
# the sqrt(4 pi)
    def _Badvection_x(field,data):
        return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'Bx') 
    add_field('Badvection_x',function=_Badvection_x, validators=std_validators, take_log=False, units=magnetic_rate_unit,sampling_type='cell')
    def _Badvection_y(field,data):
        return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'By')
    add_field('Badvection_y',function=_Badvection_y, validators=std_validators, take_log=False, units=magnetic_rate_unit,sampling_type='cell')
    def _Badvection_z(field,data):
        return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'Bz')
    add_field('Badvection_z',function=_Badvection_z, validators=std_validators, take_log=False, units=magnetic_rate_unit,sampling_type='cell')

    def _Bcompression_x(field,data):
        div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
        return -1*data['Bx']*div_v
    add_field('Bcompression_x',function=_Bcompression_x, validators=std_validators, take_log=False, units=magnetic_rate_unit,sampling_type='cell')
    def _Bcompression_y(field,data):
        div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
        return -1*data['By']*div_v
    add_field('Bcompression_y',function=_Bcompression_y, validators=std_validators, take_log=False, units=magnetic_rate_unit,sampling_type='cell')
    def _Bcompression_z(field,data):
        div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
        return -1*data['Bz']*div_v
    add_field('Bcompression_z',function=_Bcompression_z, validators=std_validators, take_log=False, units=magnetic_rate_unit,sampling_type='cell')

    def _Bstretching_x(field,data):
        return AdotDel(data,['Bx','By','Bz'], 'x-velocity')
    add_field('Bstretching_x',function=_Bstretching_x, validators=std_validators,take_log=False, units=magnetic_rate_unit,sampling_type='cell')
    def _Bstretching_y(field,data):
        return AdotDel(data,['Bx','By','Bz'], 'y-velocity')
    add_field('Bstretching_y',function=_Bstretching_y, validators=std_validators,take_log=False, units=magnetic_rate_unit,sampling_type='cell')
    def _Bstretching_z(field,data):
        return AdotDel(data,['Bx','By','Bz'], 'z-velocity')
    add_field('Bstretching_z',function=_Bstretching_z, validators=std_validators,take_log=False, units=magnetic_rate_unit,sampling_type='cell')

    """these are not physically meaningful."""
    def _normalized_LP(field,data):
        """As defined, Lpressure is the actual force, so momentum per time.  
        If we were using euler's (velocity) then Lpressure has an additional 1/rho. 
        Either way, normalize by momentum"""
        Lpressure = np.sqrt( data['Lpressure_x']**2 + data['Lpressure_y']**2 + data['Lpressure_z']**2 )
        return Lpressure/(data['density']*data['velocity_magnitude'])
    add_field('Normalized_LP',function=_normalized_LP, validators=std_validators,take_log=False, units='1/s',sampling_type='cell')
    def _normalized_LT(field,data):
        Ltension = np.sqrt( data['Ltension_x']**2 + data['Ltension_y']**2 + data['Ltension_z']**2 )
        return Ltension/(data['density']*data['velocity_magnitude'])
    add_field('Normalized_LT',function=_normalized_LT, validators=std_validators,take_log=False, units='1/s',sampling_type='cell')
    def _normalized_BA(field,data):
        """Magnetic field advection norm"""
        Badvection = np.sqrt(data['Badvection_x']**2+data['Badvection_y']**2+data['Badvection_z']**2)
        return Badvection/data['magnetic_field_strength']
    add_field('Normalized_BA',function=_normalized_BA, validators=std_validators,take_log=False,units='1/s',sampling_type='cell')
    def _normalized_BC(field,data):
        """Magnetic field advection norm"""
        Bcompression = np.sqrt(data['Bcompression_x']**2+data['Bcompression_y']**2+data['Bcompression_z']**2)
        return Bcompression/data['magnetic_field_strength']
    add_field('Normalized_BC',function=_normalized_BC, validators=std_validators,take_log=False,units='1/s',sampling_type='cell')
    def _normalized_BS(field,data):
        """Magnetic field advection norm"""
        Bstretching = np.sqrt(data['Bstretching_x']**2+data['Bstretching_y']**2+data['Bstretching_z']**2)
        return Bstretching/data['magnetic_field_strength']
    add_field('Normalized_BS',function=_normalized_BS, validators=std_validators,take_log=False,units='1/s',sampling_type='cell')

    """this is the physically meaningful quantity:
        norm_unit_velocity_x,
        norm_unit_velocity_y,
        etc...
        give \hat{A}/||A||
        The extra ||A|| is to normalize the forces.
        """
    def generate_normed_unit(quantity,magnitude,units):
        def _temp(field,data):
            return data[quantity]/(data[magnitude]*data[magnitude])

        field_name = 'norm_unit_'+quantity
        #print "created", field_name
        add_field(field_name, function = _temp, units=units, take_log=False,sampling_type='cell')

    for direction in 'xyz':
        generate_normed_unit('velocity_'+direction, 'velocity_magnitude', units='1/code_velocity')
        generate_normed_unit('momentum_'+direction, 'momentum_magnitude', units='cm**2*s/g')
        generate_normed_unit('magnetic_field_'+direction, 'magnetic_field_strength', units='1/code_magnetic')

    def generate_square(field_name, units): 
        def _temp(field,data):
            fields = ["%s_%s"%(field_name,s) for s in 'xyz']
            out = data[fields[0]]*data[fields[0]] \
                + data[fields[1]]*data[fields[1]] \
                + data[fields[2]]*data[fields[2]]
            return out
            
        new_field_name = '%s_squared'%field_name
        #print "created", new_field_name
        add_field(new_field_name, function = _temp, units = units,sampling_type='cell')

    def generate_full2(field_name, force_name, units): 
        def _temp(field,data):
            forces = ["%s_%s"%(force_name,s) for s in 'xyz']
            fields = ["%s_%s"%(field_name,s) for s in 'xyz']
            out = data[forces[0]]*data[fields[0]] \
                + data[forces[1]]*data[fields[1]] \
                + data[forces[2]]*data[fields[2]]
            return out
            
        new_field_name = '%s_f2'%force_name
        #print "created", new_field_name
        add_field(new_field_name, function = _temp, units = units,sampling_type='cell')
    def generate_increasing2(field_name, force_name): 
        def _temp(field,data):
            forces = ["%s_%s"%(force_name,s) for s in 'xyz']
            fields = ["%s_%s"%(field_name,s) for s in 'xyz']
            out = data[forces[0]]*data[fields[0]] \
                + data[forces[1]]*data[fields[1]] \
                + data[forces[2]]*data[fields[2]]
            out[ out < 0] = 0.0
            return out
        new_field_name = '%s_p2'%force_name
        #print "created", new_field_name
        add_field(new_field_name, function = _temp, units = 'code_magnetic**2/s',sampling_type='cell')
    def generate_decreasing2(field_name, force_name): 
        def _temp(field,data):
            forces = ["%s_%s"%(force_name,s) for s in 'xyz']
            fields = ["%s_%s"%(field_name,s) for s in 'xyz']
            out = data[forces[0]]*data[fields[0]] \
                + data[forces[1]]*data[fields[1]] \
                + data[forces[2]]*data[fields[2]]
            out[ out > 0] = 0.0
            return out
        new_field_name = '%s_m2'%force_name
        #print "created", new_field_name
        add_field(new_field_name, function = _temp, units = 'code_magnetic**2/s',sampling_type='cell')
    def generate_full(field_name, force_name): 
        def _temp(field,data):
            forces = ["%s_%s"%(force_name,s) for s in 'xyz']
            fields = ["norm_unit_%s_%s"%(field_name,s) for s in 'xyz']
            out = data[forces[0]]*data[fields[0]] \
                + data[forces[1]]*data[fields[1]] \
                + data[forces[2]]*data[fields[2]]
            return out
        new_field_name = '%s_f'%force_name
        #print "created", new_field_name
        add_field(new_field_name, function = _temp, units = '1/s',sampling_type='cell')
    def generate_increasing(field_name, force_name): 
        def _temp(field,data):
            forces = ["%s_%s"%(force_name,s) for s in 'xyz']
            fields = ["norm_unit_%s_%s"%(field_name,s) for s in 'xyz']
            out = data[forces[0]]*data[fields[0]] \
                + data[forces[1]]*data[fields[1]] \
                + data[forces[2]]*data[fields[2]]
            out[ out < 0] = 0.0
            return out
        new_field_name = '%s_p'%force_name
        #print "created", new_field_name
        add_field(new_field_name, function = _temp, units = '1/s',sampling_type='cell')
    def generate_decreasing(field_name, force_name): 
        def _temp(field,data):
            forces = ["%s_%s"%(force_name,s) for s in 'xyz']
            fields = ["norm_unit_%s_%s"%(field_name,s) for s in 'xyz']
            out = data[forces[0]]*data[fields[0]] \
                + data[forces[1]]*data[fields[1]] \
                + data[forces[2]]*data[fields[2]]
            out[ out > 0] = 0.0
            return out
        new_field_name = '%s_m'%force_name
        #print "created", new_field_name
        add_field(new_field_name, function = _temp, units = '1/s',sampling_type='cell')

    def generate_transverse(field_name, force_name): 
        def _temp(field,data):
            forces = [data["%s_%s"%(force_name,s)] for s in 'xyz']
            fields = [data["norm_unit_%s_%s"%(field_name,s)] for s in 'xyz']
            cross = Cross(fields,forces)
            out = np.sqrt(cross[0]**2+cross[1]**2+cross[2]**2)
            return out
        new_field_name = '%s_t'%force_name
        #print "created", new_field_name
        add_field(new_field_name, function = _temp, units = '1/s')

    generate_full2('magnetic_field','Bstretching',units='code_magnetic**2/s')
    generate_full2('magnetic_field','Bcompression',units='code_magnetic**2/s')
    generate_full2('magnetic_field','Badvection',units='code_magnetic**2/s')
    generate_increasing2('magnetic_field','Bstretching')
    generate_increasing2('magnetic_field','Bcompression')
    generate_increasing2('magnetic_field','Badvection')
    generate_decreasing2('magnetic_field','Bstretching')
    generate_decreasing2('magnetic_field','Bcompression')
    generate_decreasing2('magnetic_field','Badvection')
    generate_full('magnetic_field','Bstretching')
    generate_full('magnetic_field','Bcompression')
    generate_full('magnetic_field','Badvection')
    generate_full('momentum','Ltension') #momentum density ^2 per time
    generate_full('momentum','Lpressure')
    generate_full2('momentum','Ltension',units='g**2/(cm**4*s**3)')
    generate_full2('momentum','Lpressure',units='g**2/(cm**4*s**3)')
    generate_full2('momentum','Ltension',units='g**2/(cm**4*s**3)')
    generate_full2('momentum','Lpressure',units='g**2/(cm**4*s**3)')
    generate_square('Ltension',units='g**2/(cm**4*s**4)')
    generate_square('Lpressure',units='g**2/(cm**4*s**4)')
    generate_increasing('magnetic_field','Bstretching')
    generate_increasing('magnetic_field','Bcompression')
    generate_increasing('magnetic_field','Badvection')
    generate_increasing('momentum','Ltension')
    generate_increasing('momentum','Lpressure')
    generate_decreasing('magnetic_field','Bstretching')
    generate_decreasing('magnetic_field','Bcompression')
    generate_decreasing('magnetic_field','Badvection')
    generate_decreasing('momentum','Ltension')
    generate_decreasing('momentum','Lpressure')
    generate_transverse('magnetic_field','Bstretching')
    generate_transverse('magnetic_field','Bcompression')
    generate_transverse('magnetic_field','Badvection')
    generate_transverse('momentum','Ltension')
    generate_transverse('momentum','Lpressure')
