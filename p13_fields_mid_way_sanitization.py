mean_density = 1.0
import numpy as na
mhd = True
print "Watch units."
#def CellMassNorm(field, data):
    #return data["cell_mass"]/mean_density
#
#yt.add_field("CellMassNorm", function=CellMassNorm, display_name=r"m/m_{\rm tot}")

if 1:  #the switch                                                                                                
    def CellVolumeNorm(field, data):
        return data["cell_volume"].in_units("code_length**3")

    yt.add_field("CellVolumeNorm", function=CellVolumeNorm, display_name=r"V/L^3", units="code_length**3")

if 1:  #the switch                                                                                                

    def Delta(field, data):
        den = data['density'].in_units('code_density')
        return den/ds.arr(mean_density,"code_density")

    yt.add_field("Delta", function=Delta, display_name=r"1+\delta",units="")
if 1:  #the switch                                                                                                

    def Gravity(field, data):
        g = data.ds.arr(data.ds['GravitationalConstant'],'1/(code_density*code_time**2)')
        den = data['density'].in_units('code_density')
        md = data.ds.quan(mean_density,"code_density")
        return data['dx']*data['dx']*g*(den - md)*den

    yt.add_field("Gravity", function=Gravity, take_log=False, display_name=r"4\pi G\rho(\rho-\rho_{0})\Delta^2", units="erg/cm**3")

if 1:  #the switch                                                                                                

# ========== THERMAL SUPPORT ==========
    def KineticEnergyDensity(field, data):

        den = data['density'].in_units('code_density')
        vx = data['x-velocity'].in_units('code_length/code_time').v
        vy = data['y-velocity'].in_units('code_length/code_time').v
        vz = data['z-velocity'].in_units('code_length/code_time').v
        return 0.5*den*(data["x-velocity"]*data["x-velocity"] + \
                                    data["y-velocity"]*data["y-velocity"] + \
                                    data["z-velocity"]*data["z-velocity"])

    yt.add_field("KineticEnergyDensity", function=KineticEnergyDensity, display_name=r"\frac{1}{2}\rho v^2",units="erg/cm**3")

#   def InternalEnergy(field, data):
#       if not mhd:
#           return (data["total_energy"] - 0.5*(data["x-velocity"]*data["x-velocity"] + \
#                                               data["y-velocity"]*data["y-velocity"] + \
#                                               data["z-velocity"]*data["z-velocity"]))
#       else:
#           raise Exception("NO INTERNAL ENERGY")

#   yt.add_field("InternalEnergy", function=InternalEnergy,units=None)
if 1:  #the switch                                                                                                

    def Pressure(field, data):
        cs2 = data.ds.quan(1.0,"code_length**2/code_time**2")
        den = data["density"].in_units('code_density')
        tmp = cs2*den
        if not mhd:
            tmp *= 0.001*data["InternalEnergy"]
        return tmp

    yt.add_field("Pressure", function=Pressure, display_name=r"P",units="erg/cm**3")

if 1: #the switch

    def GradProductDensityPressure(field, data):
        # We need to set up stencils
        sl_extr_left  = slice(None,-4,None)
        sl_left       = slice(1,-3,None)
        sl_cent       = slice(2,-2,None)
        sl_right      = slice(3,-1,None)
        sl_extr_right = slice(4,None,None)
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

if 1:  #the switch                                                                                                
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

if 1:  #the switch                                                                                                
    def Lambda_therm(field, data):
        den = data['density'].in_units('code_density')
        tmp = -data["LaplacePressure"]
        tmp += data["GradProductDensityPressure"]/den
        return tmp

    yt.add_field("Lambda_therm", function=Lambda_therm, take_log=False, display_name=r"\Delta^2\!\rho\Lambda_{\rm therm}",units=None)
if 1:  #the switch                                                                                                

    def Lambda_therm_pos(field, data):
        tmp = data["Lambda_therm"]
        return na.where(tmp > 0.0, tmp, 0.0)

    yt.add_field("Lambda_therm_pos", function=Lambda_therm_pos, display_name=r"\Delta^2\!\rho\Lambda_{\rm therm\,+}",units=None)

    def Lambda_therm_neg(field, data):
        tmp = -data["Lambda_therm"]
        return na.where(tmp > 0.0, tmp, 0.0)

    yt.add_field("Lambda_therm_neg", function=Lambda_therm_neg, display_name=r"\Delta^2\!\rho\Lambda_{\rm therm\,-}",units=None)

if 1: #the switch
# ========== TURBULENT SUPPORT ==========                                                                                          

    def EnstrophyDensity(field, data):
        # We need to set up stencils
        sl_extr_left  = slice(None,-4,None)
        sl_left       = slice(1,-3,None)
        sl_cent       = slice(2,-2,None)
        sl_right      = slice(3,-1,None)
        sl_extr_right = slice(4,None,None)
        fct = 12.0
        #dx = fct * (data['dx'].flat[0])
        dx = fct
        vyx = data["y-velocity"][sl_extr_left,2:-2,2:-2]/dx
        vyx -= 8.0 * data["y-velocity"][sl_left,2:-2,2:-2]/dx
        vyx += 8.0 * data["y-velocity"][sl_right,2:-2,2:-2]/dx
        vyx -= data["y-velocity"][sl_extr_right,2:-2,2:-2]/dx
        vzx = data["z-velocity"][sl_extr_left,2:-2,2:-2]/dx
        vzx -= 8.0 * data["z-velocity"][sl_left,2:-2,2:-2]/dx
        vzx += 8.0 * data["z-velocity"][sl_right,2:-2,2:-2]/dx
        vzx -= data["z-velocity"][sl_extr_right,2:-2,2:-2]/dx
        #dy = fct * (data['dy'].flat[0])
        dy = fct
        vxy = data["x-velocity"][2:-2,sl_extr_left,2:-2]/dy
        vxy -= 8.0 * data["x-velocity"][2:-2,sl_left,2:-2]/dy
        vxy += 8.0 * data["x-velocity"][2:-2,sl_right,2:-2]/dy
        vxy -= data["x-velocity"][2:-2,sl_extr_right,2:-2]/dy
        vzy = data["z-velocity"][2:-2,sl_extr_left,2:-2]/dy
        vzy -= 8.0 * data["z-velocity"][2:-2,sl_left,2:-2]/dy
        vzy += 8.0 * data["z-velocity"][2:-2,sl_right,2:-2]/dy
        vzy -= data["z-velocity"][2:-2,sl_extr_right,2:-2]/dy
        #dz = fct * (data['dz'].flat[0])
        dz = fct
        vxz = data["x-velocity"][2:-2,2:-2,sl_extr_left]/dz
        vxz -= 8.0 * data["x-velocity"][2:-2,2:-2,sl_left]/dz
        vxz += 8.0 * data["x-velocity"][2:-2,2:-2,sl_right]/dz
        vxz -= data["x-velocity"][2:-2,2:-2,sl_extr_right]/dz
        vyz = data["y-velocity"][2:-2,2:-2,sl_extr_left]/dz
        vyz -= 8.0 * data["y-velocity"][2:-2,2:-2,sl_left]/dz
        vyz += 8.0 * data["y-velocity"][2:-2,2:-2,sl_right]/dz
        vyz -= data["y-velocity"][2:-2,2:-2,sl_extr_right]/dz
        den = data['density'].in_units('code_density')
        new_field = den
        new_field[2:-2,2:-2,2:-2] *= 0.5*((vzy-vyz)*(vzy-vyz) + (vxz-vzx)*(vxz-vzx) + (vyx-vxy)*(vyx-vxy))
        return new_field

    yt.add_field("EnstrophyDensity", function=EnstrophyDensity,
              validators=[yt.ValidateSpatial(2,["density","x-velocity","y-velocity","z-velocity"])],
              display_name=r"\frac{1}{2}\Delta^2\!\rho\omega^2",units=None)

if 1:  #the switch                                                                                                
    def RateOfStrainSqr(field, data):
        # We need to set up stencils
        sl_extr_left  = slice(None,-4,None)
        sl_left       = slice(1,-3,None)
        sl_cent       = slice(2,-2,None)
        sl_right      = slice(3,-1,None)
        sl_extr_right = slice(4,None,None)
        fct = 12.0
        dx = fct
        den = data['density'].in_units('code_density')
        vxx = data["x-velocity"][sl_extr_left,2:-2,2:-2]/dx
        vxx -= 8.0 * data["x-velocity"][sl_left,2:-2,2:-2]/dx
        vxx += 8.0 * data["x-velocity"][sl_right,2:-2,2:-2]/dx
        vxx -= data["x-velocity"][sl_extr_right,2:-2,2:-2]/dx
        vyx = data["y-velocity"][sl_extr_left,2:-2,2:-2]/dx
        vyx -= 8.0 * data["y-velocity"][sl_left,2:-2,2:-2]/dx
        vyx += 8.0 * data["y-velocity"][sl_right,2:-2,2:-2]/dx
        vyx -= data["y-velocity"][sl_extr_right,2:-2,2:-2]/dx
        vzx = data["z-velocity"][sl_extr_left,2:-2,2:-2]/dx
        vzx -= 8.0 * data["z-velocity"][sl_left,2:-2,2:-2]/dx
        vzx += 8.0 * data["z-velocity"][sl_right,2:-2,2:-2]/dx
        vzx -= data["z-velocity"][sl_extr_right,2:-2,2:-2]/dx
        dy = fct
        vxy = data["x-velocity"][2:-2,sl_extr_left,2:-2]/dy
        vxy -= 8.0 * data["x-velocity"][2:-2,sl_left,2:-2]/dy
        vxy += 8.0 * data["x-velocity"][2:-2,sl_right,2:-2]/dy
        vxy -= data["x-velocity"][2:-2,sl_extr_right,2:-2]/dy
        vyy = data["y-velocity"][2:-2,sl_extr_left,2:-2]/dy
        vyy -= 8.0 * data["y-velocity"][2:-2,sl_left,2:-2]/dy
        vyy += 8.0 * data["y-velocity"][2:-2,sl_right,2:-2]/dy
        vyy -= data["y-velocity"][2:-2,sl_extr_right,2:-2]/dy
        vzy = data["z-velocity"][2:-2,sl_extr_left,2:-2]/dy
        vzy -= 8.0 * data["z-velocity"][2:-2,sl_left,2:-2]/dy
        vzy += 8.0 * data["z-velocity"][2:-2,sl_right,2:-2]/dy
        vzy -= data["z-velocity"][2:-2,sl_extr_right,2:-2]/dy
        dz = fct
        vxz = data["x-velocity"][2:-2,2:-2,sl_extr_left]/dz
        vxz -= 8.0 * data["x-velocity"][2:-2,2:-2,sl_left]/dz
        vxz += 8.0 * data["x-velocity"][2:-2,2:-2,sl_right]/dz
        vxz -= data["x-velocity"][2:-2,2:-2,sl_extr_right]/dz
        vyz = data["y-velocity"][2:-2,2:-2,sl_extr_left]/dz
        vyz -= 8.0 * data["y-velocity"][2:-2,2:-2,sl_left]/dz
        vyz += 8.0 * data["y-velocity"][2:-2,2:-2,sl_right]/dz
        vyz -= data["y-velocity"][2:-2,2:-2,sl_extr_right]/dz
        vzz = data["z-velocity"][2:-2,2:-2,sl_extr_left]/dz
        vzz -= 8.0 * data["z-velocity"][2:-2,2:-2,sl_left]/dz
        vzz += 8.0 * data["z-velocity"][2:-2,2:-2,sl_right]/dz
        vzz -= data["z-velocity"][2:-2,2:-2,sl_extr_right]/dz
        new_field = den
        new_field[2:-2,2:-2,2:-2] *= 0.5*(2.0*(vxx*vxx + vyy*vyy + vzz*vzz) + \
                                          (vxy+vyx)*(vxy+vyx) + (vyz+vzy)*(vyz+vzy) + (vxz+vzx)*(vxz+vzx))
        return new_field

    yt.add_field("RateOfStrainSqr", function=RateOfStrainSqr,
              validators=[yt.ValidateSpatial(2,["density","x-velocity","y-velocity","z-velocity"])],
              display_name=r"\frac{1}{2}\Delta^2\!\rho|S|^2",units=None)

    def Denstrophy(field, data):
        # We need to set up stencils
        sl_extr_left  = slice(None,-4,None)
        sl_left       = slice(1,-3,None)
        sl_cent       = slice(2,-2,None)
        sl_right      = slice(3,-1,None)
        sl_extr_right = slice(4,None,None)
        fct = 12.0
        # mass-weighted velocity
        den = data['density'].in_units('code_density')
        vel_mw = [ ]
        for i in range(3):
            vel_mw.append(na.sqrt(den))
        vel_mw[0] *= data["x-velocity"]   
        vel_mw[1] *= data["y-velocity"]
        vel_mw[2] *= data["z-velocity"]
        #dx = fct * (data['dx'].flat[0])
        dx = fct
        vyx = -vel_mw[1][sl_extr_left,2:-2,2:-2]/dx
        vyx += 8.0 * vel_mw[1][sl_left,2:-2,2:-2]/dx
        vyx -= 8.0 * vel_mw[1][sl_right,2:-2,2:-2]/dx
        vyx += vel_mw[1][sl_extr_right,2:-2,2:-2]/dx
        vzx = -vel_mw[2][sl_extr_left,2:-2,2:-2]/dx
        vzx += 8.0 * vel_mw[2][sl_left,2:-2,2:-2]/dx
        vzx -= 8.0 * vel_mw[2][sl_right,2:-2,2:-2]/dx
        vzx += vel_mw[2][sl_extr_right,2:-2,2:-2]/dx
        #dy = fct * (data['dy'].flat[0])
        dy = fct
        vxy = -vel_mw[0][2:-2,sl_extr_left,2:-2]/dy
        vxy += 8.0 * vel_mw[0][2:-2,sl_left,2:-2]/dy
        vxy -= 8.0 * vel_mw[0][2:-2,sl_right,2:-2]/dy
        vxy += vel_mw[0][2:-2,sl_extr_right,2:-2]/dy
        vzy = -vel_mw[2][2:-2,sl_extr_left,2:-2]/dy
        vzy += 8.0 * vel_mw[2][2:-2,sl_left,2:-2]/dy
        vzy -= 8.0 * vel_mw[2][2:-2,sl_right,2:-2]/dy
        vzy += vel_mw[2][2:-2,sl_extr_right,2:-2]/dy
        #dz = fct * (data['dz'].flat[0])
        dz = fct
        vxz = -vel_mw[0][2:-2,2:-2,sl_extr_left]/dz
        vxz += 8.0 * vel_mw[0][2:-2,2:-2,sl_left]/dz
        vxz -= 8.0 * vel_mw[0][2:-2,2:-2,sl_right]/dz
        vxz += vel_mw[0][2:-2,2:-2,sl_extr_right]/dz
        vyz = -vel_mw[1][2:-2,2:-2,sl_extr_left]/dz
        vyz += 8.0 * vel_mw[1][2:-2,2:-2,sl_left]/dz
        vyz -= 8.0 * vel_mw[1][2:-2,2:-2,sl_right]/dz
        vyz += vel_mw[1][2:-2,2:-2,sl_extr_right]/dz
        #negative velocity gradients, but sign cancels out
        new_field = na.zeros(den.shape, dtype='float64')
        new_field[2:-2,2:-2,2:-2] = 0.5*((vzy-vyz)*(vzy-vyz) + (vxz-vzx)*(vxz-vzx) + (vyx-vxy)*(vyx-vxy))
        return new_field

    yt.add_field("Denstrophy", function=Denstrophy,
              validators=[yt.ValidateSpatial(2,["density","x-velocity","y-velocity","z-velocity"])],
              display_name=r"\Delta^2\Omega_{1/2}",units=None)
                                                                                                        
if 1:  #the switch                                                                                                
    def DiffOmegaSqrStrainSqr(field, data):
        # We need to set up stencils
        # forget about HydroMethod =2!
        sl_extr_left  = slice(None,-4,None)
        sl_left       = slice(1,-3,None)
        sl_cent       = slice(2,-2,None)
        sl_right      = slice(3,-1,None)
        sl_extr_right = slice(4,None,None)
        fct = 12.0
        dx = fct
        vxx = data["x-velocity"][sl_extr_left,2:-2,2:-2]/dx
        vxx -= 8.0 * data["x-velocity"][sl_left,2:-2,2:-2]/dx
        vxx += 8.0 * data["x-velocity"][sl_right,2:-2,2:-2]/dx
        vxx -= data["x-velocity"][sl_extr_right,2:-2,2:-2]/dx
        vyx = data["y-velocity"][sl_extr_left,2:-2,2:-2]/dx
        vyx -= 8.0 * data["y-velocity"][sl_left,2:-2,2:-2]/dx
        vyx += 8.0 * data["y-velocity"][sl_right,2:-2,2:-2]/dx
        vyx -= data["y-velocity"][sl_extr_right,2:-2,2:-2]/dx
        vzx = data["z-velocity"][sl_extr_left,2:-2,2:-2]/dx
        vzx -= 8.0 * data["z-velocity"][sl_left,2:-2,2:-2]/dx
        vzx += 8.0 * data["z-velocity"][sl_right,2:-2,2:-2]/dx
        vzx -= data["z-velocity"][sl_extr_right,2:-2,2:-2]/dx
        dy = fct
        vxy = data["x-velocity"][2:-2,sl_extr_left,2:-2]/dy
        vxy -= 8.0 * data["x-velocity"][2:-2,sl_left,2:-2]/dy
        vxy += 8.0 * data["x-velocity"][2:-2,sl_right,2:-2]/dy
        vxy -= data["x-velocity"][2:-2,sl_extr_right,2:-2]/dy
        vyy = data["y-velocity"][2:-2,sl_extr_left,2:-2]/dy
        vyy -= 8.0 * data["y-velocity"][2:-2,sl_left,2:-2]/dy
        vyy += 8.0 * data["y-velocity"][2:-2,sl_right,2:-2]/dy
        vyy -= data["y-velocity"][2:-2,sl_extr_right,2:-2]/dy
        vzy = data["z-velocity"][2:-2,sl_extr_left,2:-2]/dy
        vzy -= 8.0 * data["z-velocity"][2:-2,sl_left,2:-2]/dy
        vzy += 8.0 * data["z-velocity"][2:-2,sl_right,2:-2]/dy
        vzy -= data["z-velocity"][2:-2,sl_extr_right,2:-2]/dy
        dz = fct
        vxz = data["x-velocity"][2:-2,2:-2,sl_extr_left]/dz
        vxz -= 8.0 * data["x-velocity"][2:-2,2:-2,sl_left]/dz
        vxz += 8.0 * data["x-velocity"][2:-2,2:-2,sl_right]/dz
        vxz -= data["x-velocity"][2:-2,2:-2,sl_extr_right]/dz
        vyz = data["y-velocity"][2:-2,2:-2,sl_extr_left]/dz
        vyz -= 8.0 * data["y-velocity"][2:-2,2:-2,sl_left]/dz
        vyz += 8.0 * data["y-velocity"][2:-2,2:-2,sl_right]/dz
        vyz -= data["y-velocity"][2:-2,2:-2,sl_extr_right]/dz
        vzz = data["z-velocity"][2:-2,2:-2,sl_extr_left]/dz
        vzz -= 8.0 * data["z-velocity"][2:-2,2:-2,sl_left]/dz
        vzz += 8.0 * data["z-velocity"][2:-2,2:-2,sl_right]/dz
        vzz -= data["z-velocity"][2:-2,2:-2,sl_extr_right]/dz
        new_field = na.zeros(data["x-velocity"].shape, dtype='float64')
        # vorticiy squared
        new_field[2:-2,2:-2,2:-2] = (vzy-vyz)*(vzy-vyz) + (vxz-vzx)*(vxz-vzx) + (vyx-vxy)*(vyx-vxy)
        # strain squared
        new_field[2:-2,2:-2,2:-2] -= (2.0*(vxx*vxx + vyy*vyy + vzz*vzz) + (vxy+vyx)*(vxy+vyx) + (vyz+vzy)*(vyz+vzy) + (vxz+vzx)*(vxz+vzx))
        return new_field

    yt.add_field("DiffOmegaSqrStrainSqr", function=DiffOmegaSqrStrainSqr,
              validators=[yt.ValidateSpatial(2,["x-velocity","y-velocity","z-velocity"])],units=None)

    def Lambda_turb(field, data):
        den = data['density'].in_units('code_density')
        return 0.5*den*data["DiffOmegaSqrStrainSqr"]

    yt.add_field("Lambda_turb", function=Lambda_turb, take_log=False, display_name=r"\Delta^2\!\rho\Lambda_{\rm turb}",units=None)
if 1:  #the switch                                                                                                

    def Lambda_turb_pos(field, data):
        tmp = data["Lambda_turb"]
        return na.where(tmp > 0.0, tmp, 0.0)

    yt.add_field("Lambda_turb_pos", function=Lambda_turb_pos, display_name=r"\Delta^2\!\rho\Lambda_{\rm turb\,+}",units=None)

    def Lambda_turb_neg(field, data):
        tmp = -data["Lambda_turb"]
        return na.where(tmp > 0.0, tmp, 0.0)

    yt.add_field("Lambda_turb_neg", function=Lambda_turb_neg, display_name=r"\Delta^2\!\rho\Lambda_{\rm turb\,-}",units=None)


# ========== MAGNETIC SUPPORT ==========

    NoFourPi = 1
    def MagneticPressure(field, data):
        bx = data['Bx'].in_units('code_magnetic').v
        by = data['By'].in_units('code_magnetic').v
        bz = data['Bz'].in_units('code_magnetic').v
        if mhd:
            return (bx*bx+by*by+bz*bz)/(2*NoFourPi) #code magnetic does not get 4pi
        else:
            return na.zeros(data["Pressure"].shape)

#yt.add_field("MagneticPressure", function=MagneticPressure, display_name=r"P_{\rm magn}")
    yt.add_field("MagneticPressure", function=MagneticPressure, display_name=r"B^2\!/8\pi",units=None)

    def GradProductDensityMagneticPressure(field, data):
        # We need to set up stencils
        sl_extr_left  = slice(None,-4,None)
        sl_left       = slice(1,-3,None)
        sl_cent       = slice(2,-2,None)
        sl_right      = slice(3,-1,None)
        sl_extr_right = slice(4,None,None)
        den = data['density'].in_units('code_density')
        new_field = na.zeros(den.shape, dtype='float64')
        fct = 12.0
        #dx = fct * (data['dx'].flat[0])
        dx = fct
        f = -den[sl_extr_left,2:-2,2:-2]/dx
        f += 8.0 * den[sl_left,2:-2,2:-2]/dx
        f -= 8.0 * den[sl_right,2:-2,2:-2]/dx
        f += den[sl_extr_right,2:-2,2:-2]/dx
        g = -data["MagneticPressure"][sl_extr_left,2:-2,2:-2]/dx
        g += 8.0 * data["MagneticPressure"][sl_left,2:-2,2:-2]/dx
        g -= 8.0 * data["MagneticPressure"][sl_right,2:-2,2:-2]/dx
        g += data["MagneticPressure"][sl_extr_right,2:-2,2:-2]/dx
        new_field[2:-2,2:-2,2:-2] += f*g
        #dy = fct * (data['dy'].flat[0])
        dy = fct
        f = -den[2:-2,sl_extr_left,2:-2]/dy
        f += 8.0 * den[2:-2,sl_left,2:-2]/dy
        f -= 8.0 * den[2:-2,sl_right,2:-2]/dy
        f += den[2:-2,sl_extr_right,2:-2]/dy
        g = -data["MagneticPressure"][2:-2,sl_extr_left,2:-2]/dy
        g += 8.0 * data["MagneticPressure"][2:-2,sl_left,2:-2]/dy
        g -= 8.0 * data["MagneticPressure"][2:-2,sl_right,2:-2]/dy
        g += data["MagneticPressure"][2:-2,sl_extr_right,2:-2]/dy
        new_field[2:-2,2:-2,2:-2] += f*g
        #dz = fct * (data['dz'].flat[0])
        dz = fct
        f = -den[2:-2,2:-2,sl_extr_left]/dz
        f += 8.0 * den[2:-2,2:-2,sl_left]/dz
        f -= 8.0 * den[2:-2,2:-2,sl_right]/dz
        f += den[2:-2,2:-2,sl_extr_right]/dz
        g = -data["MagneticPressure"][2:-2,2:-2,sl_extr_left]/dz
        g += 8.0 * data["MagneticPressure"][2:-2,2:-2,sl_left]/dz
        g -= 8.0 * data["MagneticPressure"][2:-2,2:-2,sl_right]/dz
        g += data["MagneticPressure"][2:-2,2:-2,sl_extr_right]/dz
        # f and g are negative gradients, but sign cancels in the product
        new_field[2:-2,2:-2,2:-2] += f*g
        return new_field

    yt.add_field("GradProductDensityMagneticPressure", function=GradProductDensityMagneticPressure,
              validators=[yt.ValidateSpatial(2,["density","MagneticPressure"])],units=None)

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
        #ds = fct * (data['dx'].flat[0])**2.0
        f = -data["MagneticPressure"][sl_extr_left,2:-2,2:-2]/ds
        f += 16.0 * data["MagneticPressure"][sl_left,2:-2,2:-2]/ds
        f -= 30.0 * data["MagneticPressure"][sl_cent,2:-2,2:-2]/ds
        f += 16.0 * data["MagneticPressure"][sl_right,2:-2,2:-2]/ds
        f -= data["MagneticPressure"][sl_extr_right,2:-2,2:-2]/ds
        #ds = fct * (data['dy'].flat[0])**2.0
        f -= data["MagneticPressure"][2:-2,sl_extr_left,2:-2]/ds
        f += 16.0 * data["MagneticPressure"][2:-2,sl_left,2:-2]/ds
        f -= 30.0 * data["MagneticPressure"][2:-2,sl_cent,2:-2]/ds
        f += 16.0 * data["MagneticPressure"][2:-2,sl_right,2:-2]/ds
        f -= data["MagneticPressure"][2:-2,sl_extr_right,2:-2]/ds
        #ds = fct * (data['dz'].flat[0])**2.0
        f -= data["MagneticPressure"][2:-2,2:-2,sl_extr_left]/ds
        f += 16.0 * data["MagneticPressure"][2:-2,2:-2,sl_left]/ds
        f -= 30.0 * data["MagneticPressure"][2:-2,2:-2,sl_cent]/ds
        f += 16.0 * data["MagneticPressure"][2:-2,2:-2,sl_right]/ds
        f -= data["MagneticPressure"][2:-2,2:-2,sl_extr_right]/ds
        new_field = na.zeros(data["MagneticPressure"].shape, dtype='float64')
        new_field[2:-2,2:-2,2:-2] = f
        return new_field 

    yt.add_field("LaplaceMagneticPressure", function=LaplaceMagneticPressure,
              validators=[yt.ValidateSpatial(2,["MagneticPressure"])],units=None)

if 1:  #the switch                                                                                                
    def MagneticFieldCrossContractions(field, data):
        # We need to set up stencils
        sl_extr_left  = slice(None,-4,None)
        sl_left       = slice(1,-3,None)
        sl_cent       = slice(2,-2,None)
        sl_right      = slice(3,-1,None)
        sl_extr_right = slice(4,None,None)
        den = data['density'].in_units('code_density')
        new_field = na.zeros(den.shape, dtype='float64')

        # density gradient
        dens_grad = [ ]
        for i in range(3):
            dens_grad.append(na.zeros(den.shape, dtype='float64'))
        fct = 12.0
        #dx = fct * (data['dx'].flat[0])
        dx = fct
        f = -den[sl_extr_left,2:-2,2:-2]/dx
        f += 8.0 * den[sl_left,2:-2,2:-2]/dx
        f -= 8.0 * den[sl_right,2:-2,2:-2]/dx
        f += den[sl_extr_right,2:-2,2:-2]/dx
        dens_grad[0][2:-2,2:-2,2:-2] = -f
        #dy = fct * (data['dy'].flat[0])
        dy = fct
        f = -den[2:-2,sl_extr_left,2:-2]/dy
        f += 8.0 * den[2:-2,sl_left,2:-2]/dy
        f -= 8.0 * den[2:-2,sl_right,2:-2]/dy
        f += den[2:-2,sl_extr_right,2:-2]/dy
        dens_grad[1][2:-2,2:-2,2:-2] = -f
        #dz = fct * (data['dz'].flat[0])
        dz = fct
        f = -den[2:-2,2:-2,sl_extr_left]/dz
        f += 8.0 * den[2:-2,2:-2,sl_left]/dz
        f -= 8.0 * den[2:-2,2:-2,sl_right]/dz
        f += den[2:-2,2:-2,sl_extr_right]/dz
        dens_grad[2][2:-2,2:-2,2:-2] = -f

        bx = data['Bx'].in_units('code_magnetic').v
        by = data['By'].in_units('code_magnetic').v
        bz = data['Bz'].in_units('code_magnetic').v
        B = [bx,by,bz]

        Bgrad = [[ ],[ ],[ ]]
        for i in range(3):
            for j in range(3):
                Bgrad[i].append(na.zeros(den.shape, dtype='float64'))
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
                    B[j][2:-2,2:-2,2:-2]*dens_grad[i][2:-2,2:-2,2:-2]/den[2:-2,2:-2,2:-2])* \
                    Bgrad[i][j][2:-2,2:-2,2:-2]

        return new_field

    yt.add_field("MagneticFieldCrossContractions", function=MagneticFieldCrossContractions,
              validators=[yt.ValidateSpatial(2,["density","Bx","By","Bz"])],units=None)

if 1:  #the switch                                                                                                
    def Lambda_magn(field, data):
        den = data['density'].in_units('code_density')
        tmp = -data["LaplaceMagneticPressure"]
        tmp += data["GradProductDensityMagneticPressure"]/den
        tmp += data["MagneticFieldCrossContractions"]
        return tmp

    yt.add_field("Lambda_magn", function=Lambda_magn, take_log=False, display_name=r"\Delta^2\!\rho\Lambda_{\rm magn}",units=None)

if 1: #the switch

    def Lambda_magn_pos(field, data):
        tmp = data["Lambda_magn"]
        return na.where(tmp > 0.0, tmp, 0.0)

    yt.add_field("Lambda_magn_pos", function=Lambda_magn_pos, display_name=r"\Delta^2\!\rho\Lambda_{\rm magn\,+}",units=None)

    def Lambda_magn_neg(field, data):
        tmp = -data["Lambda_magn"]
        return na.where(tmp > 0.0, tmp, 0.0)

    yt.add_field("Lambda_magn_neg", function=Lambda_magn_neg, display_name=r"\Delta^2\!\rho\Lambda_{\rm magn\,-}",units=None)


# ========== TOTAL SUPPORT ==========

    def Lambda(field, data):
        if mhd:
            return (data["Lambda_therm"] + data["Lambda_turb"] + data["Lambda_magn"])
        else:
            return (data["Lambda_therm"] + data["Lambda_turb"])

#yt.add_field("Lambda", function=Lambda, display_name=r"\Delta^2\!\rho\Lambda", take_log=False,units=None)


