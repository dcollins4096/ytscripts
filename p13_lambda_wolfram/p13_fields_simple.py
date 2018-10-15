
def GradProductDensityPressure_simple(field, data):
    # We need to set up stencils
    sl_left       = slice(None,-2,None)
    sl_cent       = slice(1,-r,None)
    sl_right      = slice(2,None,None)
    den = data['density'].in_units('code_density')
    new_field = na.zeros(data["Pressure"].shape, dtype='float64')
    fct = 12.0
    #dx = fct * (data['dx'].flat[0])
    dx = fct
    f = -den[sl_extr_left,2:-2,2:-2]/dx
    f += 8.0 * den[sl_left,2:-2,2:-2]/dx
    f -= 8.0 * den[sl_right,2:-2,2:-2]/dx
    f += den[sl_extr_right,2:-2,2:-2]/dx
    g = -data["Pressure"][sl_extr_left,2:-2,2:-2]/dx
    g += 8.0 * data["Pressure"][sl_left,2:-2,2:-2]/dx
    g -= 8.0 * data["Pressure"][sl_right,2:-2,2:-2]/dx
    g += data["Pressure"][sl_extr_right,2:-2,2:-2]/dx
    new_field[2:-2,2:-2,2:-2] += f*g
    #dy = fct * (data['dy'].flat[0])
    dy = fct
    f = -den[2:-2,sl_extr_left,2:-2]/dy
    f += 8.0 * den[2:-2,sl_left,2:-2]/dy
    f -= 8.0 * den[2:-2,sl_right,2:-2]/dy
    f += den[2:-2,sl_extr_right,2:-2]/dy
    g = -data["Pressure"][2:-2,sl_extr_left,2:-2]/dy
    g += 8.0 * data["Pressure"][2:-2,sl_left,2:-2]/dy
    g -= 8.0 * data["Pressure"][2:-2,sl_right,2:-2]/dy
    g += data["Pressure"][2:-2,sl_extr_right,2:-2]/dy
    new_field[2:-2,2:-2,2:-2] += f*g
    #dz = fct * (data['dz'].flat[0])
    dz = fct
    f = -den[2:-2,2:-2,sl_extr_left]/dz
    f += 8.0 * den[2:-2,2:-2,sl_left]/dz
    f -= 8.0 * den[2:-2,2:-2,sl_right]/dz
    f += den[2:-2,2:-2,sl_extr_right]/dz
    g = -data["Pressure"][2:-2,2:-2,sl_extr_left]/dz
    g += 8.0 * data["Pressure"][2:-2,2:-2,sl_left]/dz
    g -= 8.0 * data["Pressure"][2:-2,2:-2,sl_right]/dz
    g += data["Pressure"][2:-2,2:-2,sl_extr_right]/dz
    # f and g are negative gradients, but sign cancels in the product 
    new_field[2:-2,2:-2,2:-2] += f*g
    #pdb.set_trace()
    return new_field

yt.add_field("GradProductDensityPressure", function=GradProductDensityPressure,
          validators=[yt.ValidateSpatial(2,["density","Pressure"])],units=None)

if switch:  #the switch                                                                                                
def LaplacePressure(field, data):
    # We need to set up stencils
    # forget about HydroMethod =2!
    sl_extr_left  = slice(None,-4,None)
    sl_left       = slice(1,-3,None)
    sl_cent       = slice(2,-2,None)
    sl_right      = slice(3,-1,None)
    sl_extr_right = slice(4,None,None)
    fct = 12.0
    ds = fct    
    #ds = fct * (data['dx'].flat[0])**2.0
    f  = -data["Pressure"][sl_extr_left,2:-2,2:-2]/ds
    f += 16.0 * data["Pressure"][sl_left,2:-2,2:-2]/ds
    f -= 30.0 * data["Pressure"][sl_cent,2:-2,2:-2]/ds
    f += 16.0 * data["Pressure"][sl_right,2:-2,2:-2]/ds
    f -= data["Pressure"][sl_extr_right,2:-2,2:-2]/ds
    #ds = fct * (data['dy'].flat[0])**2.0
    f -= data["Pressure"][2:-2,sl_extr_left,2:-2]/ds
    f += 16.0 * data["Pressure"][2:-2,sl_left,2:-2]/ds
    f -= 30.0 * data["Pressure"][2:-2,sl_cent,2:-2]/ds
    f += 16.0 * data["Pressure"][2:-2,sl_right,2:-2]/ds
    f -= data["Pressure"][2:-2,sl_extr_right,2:-2]/ds
    #ds = fct * (data['dz'].flat[0])**2.0
    f -= data["Pressure"][2:-2,2:-2,sl_extr_left]/ds
    f += 16.0 * data["Pressure"][2:-2,2:-2,sl_left]/ds
    f -= 30.0 * data["Pressure"][2:-2,2:-2,sl_cent]/ds
    f += 16.0 * data["Pressure"][2:-2,2:-2,sl_right]/ds
    f -= data["Pressure"][2:-2,2:-2,sl_extr_right]/ds
    new_field = na.zeros(data["Pressure"].shape, dtype='float64')
    new_field[2:-2,2:-2,2:-2] = f
    return new_field 

yt.add_field("LaplacePressure", function=LaplacePressure,
          validators=[yt.ValidateSpatial(2,["Pressure"])],units=None)

if switch:  #the switch                                                                                                
def Lambda_therm(field, data):
    den = data['density'].in_units('code_density').v
    tmp = -data["LaplacePressure"]
    tmp += data["GradProductDensityPressure"]/den
    return tmp

yt.add_field("Lambda_therm", function=Lambda_therm, take_log=False, display_name=r"\Delta^2\!\rho\Lambda_{\rm therm}",units=None)

