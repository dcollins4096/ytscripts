if 'yt' in dir():
  add_field = yt.add_field
  ValidateSpatial = yt.ValidateSpatial
  use_yt3 = True
else:
  use_yt3 = False


def Pressure(field, data):
    tmp = get_density(data)
    if not mhd:
        tmp *= 0.001*data["InternalEnergy"]
    return tmp

add_field("Pressure", function=Pressure, display_name=r"P")

def Gravity_real(field, data):
    if use_yt3:
        grav_const = data.ds['GravitationalConstant'] #4 pi G
    else:
        grav_const = data.pf['GravitationalConstant'] #4 pi G
    return grav_const*get_density(data)
add_field("Gravity_real",function=Gravity_real)

def get_density(data):
    if use_yt3:
        return data['density'].v
    else:
        return data['Density']

def get_pressure(data):
    if use_yt3:
        return data['Pressure'].v
    else:
        return data['Pressure']


def LaplacePressure_Term(field, data):
    # We need to set up stencils
    # forget about HydroMethod =2!
    sl_extr_left  = slice(None,-4,None)
    sl_left       = slice(1,-3,None)
    sl_cent       = slice(2,-2,None)
    sl_right      = slice(3,-1,None)
    sl_extr_right = slice(4,None,None)
    P = data['Pressure']
    Density = get_density(data)
    #P = get_pressure(data)
    #Density = get_density(data)
    f  = -P[sl_extr_left,2:-2,2:-2]
    f += 16.0 * P[sl_left,2:-2,2:-2]
    f -= 30.0 * P[sl_cent,2:-2,2:-2]
    f += 16.0 * P[sl_right,2:-2,2:-2]
    f -= P[sl_extr_right,2:-2,2:-2]
    f -= P[2:-2,sl_extr_left,2:-2]
    f += 16.0 * P[2:-2,sl_left,2:-2]
    f -= 30.0 * P[2:-2,sl_cent,2:-2]
    f += 16.0 * P[2:-2,sl_right,2:-2]
    f -= P[2:-2,sl_extr_right,2:-2]
    f -= P[2:-2,2:-2,sl_extr_left]
    f += 16.0 * P[2:-2,2:-2,sl_left]
    f -= 30.0 * P[2:-2,2:-2,sl_cent]
    f += 16.0 * P[2:-2,2:-2,sl_right]
    f -= P[2:-2,2:-2,sl_extr_right]
    new_field = np.zeros(P.shape, dtype='float64')
    fct = 12.0
    #ds = fct    
    ds = fct * (data['dx'].flat[0])**2.0
    new_field[2:-2,2:-2,2:-2] = f/Density[sl_cent,sl_cent,sl_cent]/ds
    return new_field 

add_field("LaplacePressure_Term", function=LaplacePressure_Term,
          validators=[ValidateSpatial(2,["Pressure","Density"])],display_name=r'$\frac{1}{\rho}\nabla^2 P$')

def GradProductDensityPressure_Term(field, data):
    # We need to set up stencils
    sl_extr_left  = slice(None,-4,None)
    sl_left       = slice(1,-3,None)
    sl_cent       = slice(2,-2,None)
    sl_right      = slice(3,-1,None)
    sl_extr_right = slice(4,None,None)
    new_field = np.zeros(data["Pressure"].shape, dtype='float64')
    fct = 12.0
    dx = fct * (data['dx'].flat[0])
    #dx = fct
    f = -data["Density"][sl_extr_left,2:-2,2:-2]/dx
    f += 8.0 * data["Density"][sl_left,2:-2,2:-2]/dx
    f -= 8.0 * data["Density"][sl_right,2:-2,2:-2]/dx
    f += data["Density"][sl_extr_right,2:-2,2:-2]/dx
    g = -data["Pressure"][sl_extr_left,2:-2,2:-2]/dx
    g += 8.0 * data["Pressure"][sl_left,2:-2,2:-2]/dx
    g -= 8.0 * data["Pressure"][sl_right,2:-2,2:-2]/dx
    g += data["Pressure"][sl_extr_right,2:-2,2:-2]/dx
    new_field[2:-2,2:-2,2:-2] += f*g
    dy = fct * (data['dy'].flat[0])
    #dy = fct
    f = -data["Density"][2:-2,sl_extr_left,2:-2]/dy
    f += 8.0 * data["Density"][2:-2,sl_left,2:-2]/dy
    f -= 8.0 * data["Density"][2:-2,sl_right,2:-2]/dy
    f += data["Density"][2:-2,sl_extr_right,2:-2]/dy
    g = -data["Pressure"][2:-2,sl_extr_left,2:-2]/dy
    g += 8.0 * data["Pressure"][2:-2,sl_left,2:-2]/dy
    g -= 8.0 * data["Pressure"][2:-2,sl_right,2:-2]/dy
    g += data["Pressure"][2:-2,sl_extr_right,2:-2]/dy
    new_field[2:-2,2:-2,2:-2] += f*g
    dz = fct * (data['dz'].flat[0])
    #dz = fct
    f = -data["Density"][2:-2,2:-2,sl_extr_left]/dz
    f += 8.0 * data["Density"][2:-2,2:-2,sl_left]/dz
    f -= 8.0 * data["Density"][2:-2,2:-2,sl_right]/dz
    f += data["Density"][2:-2,2:-2,sl_extr_right]/dz
    g = -data["Pressure"][2:-2,2:-2,sl_extr_left]/dz
    g += 8.0 * data["Pressure"][2:-2,2:-2,sl_left]/dz
    g -= 8.0 * data["Pressure"][2:-2,2:-2,sl_right]/dz
    g += data["Pressure"][2:-2,2:-2,sl_extr_right]/dz
    # f and g are negative gradients, but sign cancels in the product 
    new_field[2:-2,2:-2,2:-2] += f*g
    new_field[2:-2,2:-2,2:-2] /= (data["Density"][sl_cent,sl_cent,sl_cent]**2)
    return new_field

add_field("GradProductDensityPressure_Term", function=GradProductDensityPressure_Term,display_name =r"$\frac{1}{\rho^2}\partial_i \rho \partial_i P$",
          validators=[ValidateSpatial(2,["Density","Pressure"])])

def Lambda_therm_full(field, data):
    tmp = -data["LaplacePressure_Term"]
    tmp += data["GradProductDensityPressure_Term"]
    return tmp

add_field("Lambda_therm_full", function=Lambda_therm_full, take_log=False, display_name=r"\Lambda_{\rm therm}")

def Lambda_therm_full_per_G(field,data):
    return data['Lambda_therm_full']/data['Gravity_real']
add_field("Lambda_therm_full_per_G",function=Lambda_therm_full_per_G)


