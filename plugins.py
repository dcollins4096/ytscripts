#def _od(field,data):
#    return data['density']
#add_field('od',_od)

#def _dbg(field,data):
#    return data['DebugField']
#add_field('SFdbg',function=_dbg,take_log=False)

if 0:
    def _metal_accounting(field,data):
        output = np.zeros_like(data['density'])
        for f in ['HI_Density','HII_Density','HeI_Density','HeII_Density','HeIII_Density']:
            output += data[f]
        return output/data['density']-1.0
    add_field('metal_accounting',function=_metal_accounting,units='dimensionless',take_log=False)
    def _metal_accounting_2(field,data):
        output = np.zeros_like(data['density'])
        for f in ['HI_Density','HII_Density','HeI_Density','HeII_Density','HeIII_Density']:
            output += data[f]
        return output
    add_field('metal_accounting_2',function=_metal_accounting_2,units='g/cm**3')

    def _scaled_div_b(field,data):
        sdb = np.abs(data['enzo','DivB'])
        sdb /= data['magnetic_field_strength']
        sdb *= data.dds.max()
        return sdb
    add_field('scaled_div_b',  function=_scaled_div_b, validators=[ValidateGridType()])

    def _abs_divb(field,data):
        return np.abs(data['DivB'])
    add_field('abs_divb', function = _abs_divb)

    def _bfield_nolog(field,data):
        return data['magnetic_field_strength']
    add_field('Bstrength',function=_bfield_nolog,take_log=False,units='gauss')
    execfile('xtra_dynamo_fields.py')
    execfile('xtra_energy_fields.py')
