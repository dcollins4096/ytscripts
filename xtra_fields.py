from xtra_operators import *    
#yt.add_field('scaled_div_b',  function=_scaled_div_b, validators=[yt.ValidateGridType()])

# \dU/dt = -grad P - 1/8 pi \grad(b^2) -1/4pi B\cdot \grad 
def _od(field,data):
    return data['density'].in_units('code_density').v
yt.add_field('od',_od)
if 1:
    momentum_units = 'g/(cm**2*s)'
    def _momentum_x(field,data):
        return data['density']*data['velocity_x']
    yt.add_field('momentum_x',function=_momentum_x,units=momentum_units)
    def _momentum_y(field,data):
        return data['density']*data['velocity_y']
    yt.add_field('momentum_y',function=_momentum_y,units=momentum_units)
    def _momentum_z(field,data):
        return data['density']*data['velocity_z']
    yt.add_field('momentum_z',function=_momentum_z,units=momentum_units)
    def _momentum_magnitude(field,data):
        out = data['momentum_x']**2+data['momentum_y']**2 + data['momentum_z']**2
        return np.sqrt(out)
    yt.add_field('momentum_magnitude',function=_momentum_magnitude,units=momentum_units)

    def _mean_square_velocity(field,data):
        out = data['velocity_x']**2+data['velocity_y']**2 + data['velocity_z']**2
        return out
    yt.add_field('mean_square_velocity',function=_mean_square_velocity,units="cm**2/s**2")

    eng_units = 'g/(cm*s**2)'
    def _eng_x(field,data):
        return data['density']*data['velocity_x']*data['velocity_x']
    yt.add_field('eng_x',function=_eng_x,units=eng_units)
    def _eng_y(field,data):
        return data['density']*data['velocity_y']*data['velocity_y']
    yt.add_field('eng_y',function=_eng_y,units=eng_units)
    def _eng_z(field,data):
        return data['density']*data['velocity_z']*data['velocity_z']
    yt.add_field('eng_z',function=_eng_z,units=eng_units)

if 0:
    def _metal_accounting(field,data):
        output = np.zeros_like(data['density'])
        for f in ['HI_Density','HII_Density','HeI_Density','HeII_Density','HeIII_Density']:
            output += data[f]
        return output/data['density']-1.0
    yt.add_field('metal_accounting',function=_metal_accounting,units='dimensionless',take_log=False)
    def _metal_accounting_2(field,data):
        output = np.zeros_like(data['density'])
        for f in ['HI_Density','HII_Density','HeI_Density','HeII_Density','HeIII_Density']:
            output += data[f]
        return output
    yt.add_field('metal_accounting_2',function=_metal_accounting_2,units='g/cm**3')
if 0:
    ef('xtra_operators.py')
    def _scaled_div_b(field,data):
        sdb = np.abs(data['enzo','DivB'])
        sdb /= data['magnetic_field_strength']
        sdb *= data.dds.max()
        return sdb
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
        return na.sqrt( dBdT[0]**2+dBdT[1]**2 +dBdT[2]**2)
    yt.add_field('InductionNorm',function=_InductionNorm,
                    validators=[yt.ValidateSpatial(1,['Bx','By','Bz', 'x-velocity','y-velocity','z-velocity'])])

    std_validators = [yt.ValidateSpatial(1,['Bx','By','Bz', 'x-velocity','y-velocity','z-velocity'])]
    std_validators_2 = [yt.ValidateSpatial(1,['magnetic_field_x','magnetic_field_y','magnetic_field_z', 'x-velocity','y-velocity','z-velocity'])]
    momentum_rate_units = 'g/(cm**2*s**2)' #momentum density per time = (gram cm/s)/cm^3/time = g /(cm^2
    def _Ltension_x(field,data):
        return  AdotDel(data,['magnetic_field_x','magnetic_field_y','magnetic_field_z'],'magnetic_field_x')/(4*np.pi)
    yt.add_field('Ltension_x', function=_Ltension_x,validators=std_validators, take_log=False, units = momentum_rate_units) 
    def _Ltension_y(field,data):
        return  AdotDel(data,['magnetic_field_x','magnetic_field_y','magnetic_field_z'],'magnetic_field_y')/(4*np.pi)
    yt.add_field('Ltension_y', function=_Ltension_y,validators=std_validators, take_log=False, units = momentum_rate_units) 
    def _Ltension_z(field,data):
        return  AdotDel(data,['magnetic_field_x','magnetic_field_y','magnetic_field_z'],'magnetic_field_z')/(4*np.pi)
    yt.add_field('Ltension_z', function=_Ltension_z,validators=std_validators, take_log=False, units = momentum_rate_units) 
    def _Lpressure_x(field,data):
        return -1*grad(data,'magnetic_energy',0) #the 1/8pi is in the magnetic energy.
    yt.add_field('Lpressure_x', function=_Lpressure_x,validators=std_validators, take_log=False, units=momentum_rate_units)
    def _Lpressure_y(field,data):
        return -1*grad(data,'magnetic_energy',1)
    yt.add_field('Lpressure_y', function=_Lpressure_y,validators=std_validators, take_log=False, units=momentum_rate_units)
    def _Lpressure_z(field,data):
        return -1*grad(data,'magnetic_energy',2)
    yt.add_field('Lpressure_z', function=_Lpressure_z,validators=std_validators, take_log=False, units=momentum_rate_units)

    magnetic_rate_units = 'code_magnetic/s'
#On these units:  Using Bx and magnetic_field_x are interchangeable; the units ('code_magnetic' vs 'gauss') will take care of 
# the sqrt(4 pi)
    def _Badvection_x(field,data):
        return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'Bx') 
    yt.add_field('Badvection_x',function=_Badvection_x, validators=std_validators, take_log=False, units=magnetic_rate_units)
    def _Badvection_y(field,data):
        return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'By')
    yt.add_field('Badvection_y',function=_Badvection_y, validators=std_validators, take_log=False, units=magnetic_rate_units)
    def _Badvection_z(field,data):
        return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'Bz')
    yt.add_field('Badvection_z',function=_Badvection_z, validators=std_validators, take_log=False, units=magnetic_rate_units)

    def _Bcompression_x(field,data):
        div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
        return -1*data['Bx']*div_v
    yt.add_field('Bcompression_x',function=_Bcompression_x, validators=std_validators, take_log=False, units=magnetic_rate_units)
    def _Bcompression_y(field,data):
        div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
        return -1*data['By']*div_v
    yt.add_field('Bcompression_y',function=_Bcompression_y, validators=std_validators, take_log=False, units=magnetic_rate_units)
    def _Bcompression_z(field,data):
        div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
        return -1*data['Bz']*div_v
    yt.add_field('Bcompression_z',function=_Bcompression_z, validators=std_validators, take_log=False, units=magnetic_rate_units)

    def _Bstretching_x(field,data):
        return AdotDel(data,['Bx','By','Bz'], 'x-velocity')
    yt.add_field('Bstretching_x',function=_Bstretching_x, validators=std_validators,take_log=False, units=magnetic_rate_units)
    def _Bstretching_y(field,data):
        return AdotDel(data,['Bx','By','Bz'], 'y-velocity')
    yt.add_field('Bstretching_y',function=_Bstretching_y, validators=std_validators,take_log=False, units=magnetic_rate_units)
    def _Bstretching_z(field,data):
        return AdotDel(data,['Bx','By','Bz'], 'z-velocity')
    yt.add_field('Bstretching_z',function=_Bstretching_z, validators=std_validators,take_log=False, units=magnetic_rate_units)

    """these are not physically meaningful."""
    def _normalized_LP(field,data):
        """As defined, Lpressure is the actual force, so momentum per time.  
        If we were using euler's (velocity) then Lpressure has an additional 1/rho. 
        Either way, normalize by momentum"""
        Lpressure = np.sqrt( data['Lpressure_x']**2 + data['Lpressure_y']**2 + data['Lpressure_z']**2 )
        return Lpressure/(data['density']*data['velocity_magnitude'])
    yt.add_field('Normalized_LP',function=_normalized_LP, validators=std_validators,take_log=False, units='1/s')
    def _normalized_LT(field,data):
        Ltension = np.sqrt( data['Ltension_x']**2 + data['Ltension_y']**2 + data['Ltension_z']**2 )
        return Ltension/(data['density']*data['velocity_magnitude'])
    yt.add_field('Normalized_LT',function=_normalized_LT, validators=std_validators,take_log=False, units='1/s')
    def _normalized_BA(field,data):
        """Magnetic field advection norm"""
        Badvection = np.sqrt(data['Badvection_x']**2+data['Badvection_y']**2+data['Badvection_z']**2)
        return Badvection/data['magnetic_field_strength']
    yt.add_field('Normalized_BA',function=_normalized_BA, validators=std_validators,take_log=False,units='1/s')
    def _normalized_BC(field,data):
        """Magnetic field advection norm"""
        Bcompression = np.sqrt(data['Bcompression_x']**2+data['Bcompression_y']**2+data['Bcompression_z']**2)
        return Bcompression/data['magnetic_field_strength']
    yt.add_field('Normalized_BC',function=_normalized_BC, validators=std_validators,take_log=False,units='1/s')
    def _normalized_BS(field,data):
        """Magnetic field advection norm"""
        Bstretching = np.sqrt(data['Bstretching_x']**2+data['Bstretching_y']**2+data['Bstretching_z']**2)
        return Bstretching/data['magnetic_field_strength']
    yt.add_field('Normalized_BS',function=_normalized_BS, validators=std_validators,take_log=False,units='1/s')

    if 1:
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
            print "created", field_name
            yt.add_field(field_name, function = _temp, units=units, take_log=False)

        for direction in 'xyz':
            generate_normed_unit('velocity_'+direction, 'velocity_magnitude', units='1/code_velocity')
            generate_normed_unit('momentum_'+direction, 'momentum_magnitude', units='cm**2*s/g')
            generate_normed_unit('magnetic_field_'+direction, 'magnetic_field_strength', units='1/code_magnetic')

        def generate_full(field_name, force_name): 
            def _temp(field,data):
                forces = ["%s_%s"%(force_name,s) for s in 'xyz']
                fields = ["norm_unit_%s_%s"%(field_name,s) for s in 'xyz']
                out = data[forces[0]]*data[fields[0]] \
                    + data[forces[1]]*data[fields[1]] \
                    + data[forces[2]]*data[fields[2]]
                return out
            new_field_name = '%s_f'%force_name
            print "created", new_field_name
            yt.add_field(new_field_name, function = _temp, units = '1/s')
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
            print "created", new_field_name
            yt.add_field(new_field_name, function = _temp, units = '1/s')
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
            print "created", new_field_name
            yt.add_field(new_field_name, function = _temp, units = '1/s')
        def generate_transverse(field_name, force_name): 
            def _temp(field,data):
                forces = [data["%s_%s"%(force_name,s)] for s in 'xyz']
                fields = [data["norm_unit_%s_%s"%(field_name,s)] for s in 'xyz']
                cross = Cross(fields,forces)
                out = np.sqrt(cross[0]**2+cross[1]**2+cross[2]**2)
                return out
            new_field_name = '%s_t'%force_name
            print "created", new_field_name
            yt.add_field(new_field_name, function = _temp, units = '1/s')

        generate_full('magnetic_field','Bstretching')
        generate_full('magnetic_field','Bcompression')
        generate_full('magnetic_field','Badvection')
        generate_full('momentum','Ltension')
        generate_full('momentum','Lpressure')
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
        all_x=['Bstretching',
        'Bcompression',
        'Badvection',
        'Ltension',
        'Lpressure']



    print "MONKEY ON THE CAR"

    if 0:
        def _scaled_div_b(field,data):
            sdb = np.abs(data['enzo','DivB'])
            sdb /= data['magnetic_field_strength']
            sdb *= data.dds.max()
            return sdb
        yt.add_field('scaled_div_b',  function=_scaled_div_b, validators=[yt.ValidateGridType()])

        def _abs_divb(field,data):
            return np.abs(data['DivB'])
        yt.add_field('abs_divb', function = _abs_divb)

#print "MONKEY ON THE CAR"
