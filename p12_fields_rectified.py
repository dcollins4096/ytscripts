if 'yt' in dir():
  add_field = yt.add_field
  ValidateSpatial = yt.ValidateSpatial
  use_yt3 = True
else:
  use_yt3 = False



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

def get_N_velocity(data,N):
    fieldname = "%s-velocity"%N
    if use_yt3:
        return data[fieldname].v
    else:
        return data[fieldname]

def get_BN(data,N):
    fieldname = "B"+N
    if use_yt3:
        return data[fieldname].v*np.sqrt(4*np.pi)
    else:
        return data[fieldname]

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

    
def DiffOmegaSqrStrainSqr_d2(field, data):
    # We need to set up stencils
    # forget about HydroMethod =2!
    sl_extr_left  = slice(None,-4,None)
    sl_left       = slice(1,-3,None)
    sl_cent       = slice(2,-2,None)
    sl_right      = slice(3,-1,None)
    sl_extr_right = slice(4,None,None)
    fct = 12.0
    #dx = fct
    dx = fct * (data['dx'].flat[0])
    VX = get_N_velocity(data,"x")
    VY = get_N_velocity(data,"y")
    VZ = get_N_velocity(data,"z")
    vxx = VX[sl_extr_left,2:-2,2:-2]/dx
    vxx -= 8.0 * VX[sl_left,2:-2,2:-2]/dx
    vxx += 8.0 * VX[sl_right,2:-2,2:-2]/dx
    vxx -= VX[sl_extr_right,2:-2,2:-2]/dx
    vyx = VY[sl_extr_left,2:-2,2:-2]/dx
    vyx -= 8.0 * VY[sl_left,2:-2,2:-2]/dx
    vyx += 8.0 * VY[sl_right,2:-2,2:-2]/dx
    vyx -= VY[sl_extr_right,2:-2,2:-2]/dx
    vzx = VZ[sl_extr_left,2:-2,2:-2]/dx
    vzx -= 8.0 * VZ[sl_left,2:-2,2:-2]/dx
    vzx += 8.0 * VZ[sl_right,2:-2,2:-2]/dx
    vzx -= VZ[sl_extr_right,2:-2,2:-2]/dx
    #dy = fct
    dy = fct * (data['dy'].flat[0])
    vxy = VX[2:-2,sl_extr_left,2:-2]/dy
    vxy -= 8.0 * VX[2:-2,sl_left,2:-2]/dy
    vxy += 8.0 * VX[2:-2,sl_right,2:-2]/dy
    vxy -= VX[2:-2,sl_extr_right,2:-2]/dy
    vyy = VY[2:-2,sl_extr_left,2:-2]/dy
    vyy -= 8.0 * VY[2:-2,sl_left,2:-2]/dy
    vyy += 8.0 * VY[2:-2,sl_right,2:-2]/dy
    vyy -= VY[2:-2,sl_extr_right,2:-2]/dy
    vzy = VZ[2:-2,sl_extr_left,2:-2]/dy
    vzy -= 8.0 * VZ[2:-2,sl_left,2:-2]/dy
    vzy += 8.0 * VZ[2:-2,sl_right,2:-2]/dy
    vzy -= VZ[2:-2,sl_extr_right,2:-2]/dy
    dz = fct
    dz = fct * (data['dz'].flat[0])
    vxz = VX[2:-2,2:-2,sl_extr_left]/dz
    vxz -= 8.0 * VX[2:-2,2:-2,sl_left]/dz
    vxz += 8.0 * VX[2:-2,2:-2,sl_right]/dz
    vxz -= VX[2:-2,2:-2,sl_extr_right]/dz
    vyz = VY[2:-2,2:-2,sl_extr_left]/dz
    vyz -= 8.0 * VY[2:-2,2:-2,sl_left]/dz
    vyz += 8.0 * VY[2:-2,2:-2,sl_right]/dz
    vyz -= VY[2:-2,2:-2,sl_extr_right]/dz
    vzz = VZ[2:-2,2:-2,sl_extr_left]/dz
    vzz -= 8.0 * VZ[2:-2,2:-2,sl_left]/dz
    vzz += 8.0 * VZ[2:-2,2:-2,sl_right]/dz
    vzz -= VZ[2:-2,2:-2,sl_extr_right]/dz
    new_field = np.zeros(VX.shape, dtype='float64')
    # vorticiy squared
    new_field[2:-2,2:-2,2:-2] = (vzy-vyz)*(vzy-vyz) + (vxz-vzx)*(vxz-vzx) + (vyx-vxy)*(vyx-vxy)
    # strain squared
    new_field[2:-2,2:-2,2:-2] -= (2.0*(vxx*vxx + vyy*vyy + vzz*vzz) + (vxy+vyx)*(vxy+vyx) + (vyz+vzy)*(vyz+vzy) + (vxz+vzx)*(vxz+vzx))
    return new_field

add_field("DiffOmegaSqrStrainSqr_d2", function=DiffOmegaSqrStrainSqr_d2,
          validators=[ValidateSpatial(2,["x-velocity","y-velocity","z-velocity"])])

def Lambda_turb(field, data):
    return 0.5*data["DiffOmegaSqrStrainSqr_d2"]

add_field("Lambda_turb", function=Lambda_turb, take_log=False, display_name=r"\Lambda_{\rm turb}")

def MagneticPressure(field, data):
    if mhd:
        return (get_BN(data,"x")**2 + 
                get_BN(data,"y")**2 + 
                get_BN(data,"z")**2)/(8*math.pi)
    else:
        return np.zeros(data["Pressure"].shape)

#add_field("MagneticPressure", function=MagneticPressure, display_name=r"P_{\rm magn}")
add_field("MagneticPressure", function=MagneticPressure, display_name=r"B^2\!/8\pi")

def LaplaceMagneticPressure(field, data):
    # We need to set up stencils
    # forget about HydroMethod =2!
    sl_extr_left  = slice(None,-4,None)
    sl_left       = slice(1,-3,None)
    sl_cent       = slice(2,-2,None)
    sl_right      = slice(3,-1,None)
    sl_extr_right = slice(4,None,None)
    fct = 12.0
    ds = fct    
    ds = fct * (data['dx'].flat[0])**2.0
    f = -data["MagneticPressure"][sl_extr_left,2:-2,2:-2]/ds
    f += 16.0 * data["MagneticPressure"][sl_left,2:-2,2:-2]/ds
    f -= 30.0 * data["MagneticPressure"][sl_cent,2:-2,2:-2]/ds
    f += 16.0 * data["MagneticPressure"][sl_right,2:-2,2:-2]/ds
    f -= data["MagneticPressure"][sl_extr_right,2:-2,2:-2]/ds
    ds = fct * (data['dy'].flat[0])**2.0
    f -= data["MagneticPressure"][2:-2,sl_extr_left,2:-2]/ds
    f += 16.0 * data["MagneticPressure"][2:-2,sl_left,2:-2]/ds
    f -= 30.0 * data["MagneticPressure"][2:-2,sl_cent,2:-2]/ds
    f += 16.0 * data["MagneticPressure"][2:-2,sl_right,2:-2]/ds
    f -= data["MagneticPressure"][2:-2,sl_extr_right,2:-2]/ds
    ds = fct * (data['dz'].flat[0])**2.0
    f -= data["MagneticPressure"][2:-2,2:-2,sl_extr_left]/ds
    f += 16.0 * data["MagneticPressure"][2:-2,2:-2,sl_left]/ds
    f -= 30.0 * data["MagneticPressure"][2:-2,2:-2,sl_cent]/ds
    f += 16.0 * data["MagneticPressure"][2:-2,2:-2,sl_right]/ds
    f -= data["MagneticPressure"][2:-2,2:-2,sl_extr_right]/ds
    new_field = np.zeros(data["MagneticPressure"].shape, dtype='float64')
    new_field[2:-2,2:-2,2:-2] = f
    return new_field 

add_field("LaplaceMagneticPressure", function=LaplaceMagneticPressure,
          validators=[ValidateSpatial(2,["MagneticPressure"])])

def MagneticFieldCrossContractions(field, data):
    # We need to set up stencils
    sl_extr_left  = slice(None,-4,None)
    sl_left       = slice(1,-3,None)
    sl_cent       = slice(2,-2,None)
    sl_right      = slice(3,-1,None)
    sl_extr_right = slice(4,None,None)
    Density = get_density(data)
    new_field = np.zeros(Density.shape, dtype='float64')
    # density gradient
    dens_grad = [ ]
    for i in range(3):
        dens_grad.append(np.zeros(Density.shape, dtype='float64'))
    fct = 12.0
    dx = fct * (data['dx'].flat[0])
    #dx = fct
    Bx = get_BN(data,'x')
    By = get_BN(data,'y')
    Bz = get_BN(data,'z')
    f = -Density[sl_extr_left,2:-2,2:-2]/dx
    f += 8.0 * Density[sl_left,2:-2,2:-2]/dx
    f -= 8.0 * Density[sl_right,2:-2,2:-2]/dx
    f += Density[sl_extr_right,2:-2,2:-2]/dx
    dens_grad[0][2:-2,2:-2,2:-2] = -f
    dy = fct * (data['dy'].flat[0])
    #dy = fct
    f = -Density[2:-2,sl_extr_left,2:-2]/dy
    f += 8.0 * Density[2:-2,sl_left,2:-2]/dy
    f -= 8.0 * Density[2:-2,sl_right,2:-2]/dy
    f += Density[2:-2,sl_extr_right,2:-2]/dy
    dens_grad[1][2:-2,2:-2,2:-2] = -f
    dz = fct * (data['dz'].flat[0])
    #dz = fct
    f = -Density[2:-2,2:-2,sl_extr_left]/dz
    f += 8.0 * Density[2:-2,2:-2,sl_left]/dz
    f -= 8.0 * Density[2:-2,2:-2,sl_right]/dz
    f += Density[2:-2,2:-2,sl_extr_right]/dz
    dens_grad[2][2:-2,2:-2,2:-2] = -f

    #B = [data["Bx"],data["By"],data["Bz"]]
    B = [Bx,By,Bz]

    Bgrad = [[ ],[ ],[ ]]
    for i in range(3):
        for j in range(3):
            Bgrad[i].append(np.zeros(Density.shape, dtype='float64'))
    for i in range(3):
        # B_i,x
        f = -B[i][sl_extr_left,2:-2,2:-2]/dx
        f += 8.0 * B[i][sl_left,2:-2,2:-2]/dx
        f -= 8.0 * B[i][sl_right,2:-2,2:-2]/dx
        f += B[i][sl_extr_right,2:-2,2:-2]/dx
        Bgrad[i][0][2:-2,2:-2,2:-2] = -f
        # B_i,y
        f = -B[i][2:-2,sl_extr_left,2:-2]/dy
        f += 8.0 * B[i][2:-2,sl_left,2:-2]/dy
        f -= 8.0 * B[i][2:-2,sl_right,2:-2]/dy
        f += B[i][2:-2,sl_extr_right,2:-2]/dy
        Bgrad[i][1][2:-2,2:-2,2:-2] = -f
        # B_i,z
        f = -B[i][2:-2,2:-2,sl_extr_left]/dz
        f += 8.0 * B[i][2:-2,2:-2,sl_left]/dz
        f -= 8.0 * B[i][2:-2,2:-2,sl_right]/dz
        f += B[i][2:-2,2:-2,sl_extr_right]/dz
        Bgrad[i][2][2:-2,2:-2,2:-2] = -f

    for i in range(3):
        for j in range(3):
            new_field[2:-2,2:-2,2:-2] += (Bgrad[j][i][2:-2,2:-2,2:-2] - \
                B[j][2:-2,2:-2,2:-2]*dens_grad[i][2:-2,2:-2,2:-2]/Density[2:-2,2:-2,2:-2])* \
                Bgrad[i][j][2:-2,2:-2,2:-2]

    return new_field/(4*math.pi)

add_field("MagneticFieldCrossContractions", function=MagneticFieldCrossContractions,
          validators=[ValidateSpatial(2,["Density","Bx","By","Bz"])])

def GradProductDensityMagneticPressure(field, data):
    # We need to set up stencils
    sl_extr_left  = slice(None,-4,None)
    sl_left       = slice(1,-3,None)
    sl_cent       = slice(2,-2,None)
    sl_right      = slice(3,-1,None)
    sl_extr_right = slice(4,None,None)
    fct = 12.0
    dx = fct * (data['dx'].flat[0])
    Density = get_density(data)
    new_field = np.zeros(Density.shape, dtype='float64')
    #dx = fct
    f = -Density[sl_extr_left,2:-2,2:-2]/dx
    f += 8.0 * Density[sl_left,2:-2,2:-2]/dx
    f -= 8.0 * Density[sl_right,2:-2,2:-2]/dx
    f += Density[sl_extr_right,2:-2,2:-2]/dx
    g = -data["MagneticPressure"][sl_extr_left,2:-2,2:-2]/dx
    g += 8.0 * data["MagneticPressure"][sl_left,2:-2,2:-2]/dx
    g -= 8.0 * data["MagneticPressure"][sl_right,2:-2,2:-2]/dx
    g += data["MagneticPressure"][sl_extr_right,2:-2,2:-2]/dx
    new_field[2:-2,2:-2,2:-2] += f*g
    dy = fct * (data['dy'].flat[0])
    #dy = fct
    f = -Density[2:-2,sl_extr_left,2:-2]/dy
    f += 8.0 * Density[2:-2,sl_left,2:-2]/dy
    f -= 8.0 * Density[2:-2,sl_right,2:-2]/dy
    f += Density[2:-2,sl_extr_right,2:-2]/dy
    g = -data["MagneticPressure"][2:-2,sl_extr_left,2:-2]/dy
    g += 8.0 * data["MagneticPressure"][2:-2,sl_left,2:-2]/dy
    g -= 8.0 * data["MagneticPressure"][2:-2,sl_right,2:-2]/dy
    g += data["MagneticPressure"][2:-2,sl_extr_right,2:-2]/dy
    new_field[2:-2,2:-2,2:-2] += f*g
    dz = fct * (data['dz'].flat[0])
    #dz = fct
    f = -Density[2:-2,2:-2,sl_extr_left]/dz
    f += 8.0 * Density[2:-2,2:-2,sl_left]/dz
    f -= 8.0 * Density[2:-2,2:-2,sl_right]/dz
    f += Density[2:-2,2:-2,sl_extr_right]/dz
    g = -data["MagneticPressure"][2:-2,2:-2,sl_extr_left]/dz
    g += 8.0 * data["MagneticPressure"][2:-2,2:-2,sl_left]/dz
    g -= 8.0 * data["MagneticPressure"][2:-2,2:-2,sl_right]/dz
    g += data["MagneticPressure"][2:-2,2:-2,sl_extr_right]/dz
    # f and g are negative gradients, but sign cancels in the product
    new_field[2:-2,2:-2,2:-2] += f*g
    return new_field

add_field("GradProductDensityMagneticPressure", function=GradProductDensityMagneticPressure,
          validators=[ValidateSpatial(2,["Density","MagneticPressure"])])

def Lambda_magn(field, data):
    Density = get_density(data)
    tmp = -data["LaplaceMagneticPressure"]/Density
    #tmp = -data["LaplacianMagneticPressure_27pt"]/Density
    tmp += data["GradProductDensityMagneticPressure"]/Density**2
    tmp += data["MagneticFieldCrossContractions"]/Density
    return tmp
add_field("Lambda_magn", function=Lambda_magn, take_log=False, display_name=r"\Lambda_{\rm magn}")

def Lambda_magn_27(field, data):
    Density = get_density(data)
    tmp = -data["LaplacianMagneticPressure_27pt"]/Density
    tmp += data["GradProductDensityMagneticPressure"]/Density**2
    tmp += data["MagneticFieldCrossContractions"]/Density
    return tmp
add_field("Lambda_magn_27_bh", function=Lambda_magn_27, take_log=False, display_name=r"\Lambda_{\rm magn}")



def Lambda(field, data):
    if mhd:
        return (data["Lambda_therm_full"] + data["Lambda_turb"] + data["Lambda_magn"])
    else:
        return (data["Lambda_therm"] + data["Lambda_turb"])

add_field("Lambda", function=Lambda, display_name=r"\Lambda", take_log=False,units=None)

def make_field_per_g(field_name):
    def func(field,data):
        return data[field_name]/data['Gravity_real']
    outname = field_name+"_per_G"
    add_field(outname,function=func)

def make_field_pos(field_name):
    def func(field,data):
        tmp = data[field_name]
        return np.where(tmp > 0.0, tmp, 0.0)
    outname = field_name+"_pos"
    add_field(outname, function = func)

def make_field_neg(field_name):
    def func(field,data):
        tmp = -data[field_name]
        return np.where(tmp > 0.0, tmp, 0.0)
    outname = field_name+"_neg"
    add_field(outname, function = func)


def color_func(field_name):
    c='k'
    if field_name.startswith("Lambda_therm"):
        c='b'
    if field_name.startswith("Lambda_turb"):
        c='r'
    if field_name.startswith("Lambda_magn"):
        c='g'
    return c

def line_func(field_name):
    linestyle = '-'
    try:
        if field_name.index("neg"):
            linestyle = "--"
    except:
        linestyle = '-'
    return linestyle



make_field_per_g("Lambda_therm_full")
make_field_pos("Lambda_therm_full_per_G")
make_field_neg("Lambda_therm_full_per_G")
make_field_per_g("Lambda_turb")
make_field_pos("Lambda_turb_per_G")
make_field_neg("Lambda_turb_per_G")
make_field_per_g("Lambda_magn")
make_field_pos("Lambda_magn_per_G")
make_field_neg("Lambda_magn_per_G")


make_field_pos("Lambda_magn")
make_field_neg("Lambda_magn")
make_field_pos("Lambda_magn_27_simple")
make_field_neg("Lambda_magn_27_simple")
