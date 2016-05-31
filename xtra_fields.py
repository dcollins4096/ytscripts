from xtra_operators import *    
ef('xtra_operators.py')
def _scaled_div_b(field,data):
    sdb = np.abs(data['enzo','DivB'])
    sdb /= data['magnetic_field_strength']
    sdb *= data.dds.max()
    return sdb
#yt.add_field('scaled_div_b',  function=_scaled_div_b, validators=[yt.ValidateGridType()])

# \dU/dt = -grad P - 1/8 pi \grad(b^2) -1/4pi B\cdot \grad 

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

def _Ltension_x(field,data):
    return  AdotDel(data,['Bx','By','Bz'],'Bx')
yt.add_field('Ltension_x', function=_Ltension_x,validators=std_validators, take_log=False)
def _Ltension_y(field,data):
    return  AdotDel(data,['Bx','By','Bz'],'By')
yt.add_field('Ltension_y', function=_Ltension_y,validators=std_validators, take_log=False)
def _Ltension_z(field,data):
    return  AdotDel(data,['Bx','By','Bz'],'Bz')
yt.add_field('Ltension_z', function=_Ltension_z,validators=std_validators, take_log=False)

def _Lpressure_x(field,data):
    return -1*grad(data,'magnetic_energy',0)
yt.add_field('Lpressure_x', function=_Lpressure_x,validators=std_validators, take_log=False)
def _Lpressure_y(field,data):
    return -1*grad(data,'magnetic_energy',1)
yt.add_field('Lpressure_y', function=_Lpressure_y,validators=std_validators, take_log=False)
def _Lpressure_z(field,data):
    return -1*grad(data,'magnetic_energy',2)
yt.add_field('Lpressure_z', function=_Lpressure_z,validators=std_validators, take_log=False)


def _Badvection_x(field,data):
    return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'Bx')
yt.add_field('Badvection_x',function=_Badvection_x, validators=std_validators, take_log=False)
def _Badvection_y(field,data):
    return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'By')
yt.add_field('Badvection_y',function=_Badvection_y, validators=std_validators, take_log=False)
def _Badvection_z(field,data):
    return -1*AdotDel(data, ['x-velocity','y-velocity','z-velocity'], 'Bz')
yt.add_field('Badvection_z',function=_Badvection_z, validators=std_validators, take_log=False)

def _Bcompression_x(field,data):
    div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
    return -1*data['Bx'].v*div_v
yt.add_field('Bcompression_x',function=_Bcompression_x, validators=std_validators, take_log=False)
def _Bcompression_y(field,data):
    div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
    return -1*data['By'].v*div_v
yt.add_field('Bcompression_y',function=_Bcompression_y, validators=std_validators, take_log=False)
def _Bcompression_z(field,data):
    div_v = grad(data,'x-velocity',0)+grad(data,'y-velocity',1)+grad(data,'z-velocity',2)
    return -1*data['Bz'].v*div_v
yt.add_field('Bcompression_z',function=_Bcompression_z, validators=std_validators, take_log=False)

def _Bstretching_x(field,data):
    return AdotDel(data,['Bx','By','Bz'], 'x-velocity')
yt.add_field('Bstretching_x',function=_Bstretching_x, validators=std_validators,take_log=False)
def _Bstretching_y(field,data):
    return AdotDel(data,['Bx','By','Bz'], 'y-velocity')
yt.add_field('Bstretching_y',function=_Bstretching_y, validators=std_validators,take_log=False)
def _Bstretching_z(field,data):
    return AdotDel(data,['Bx','By','Bz'], 'z-velocity')
yt.add_field('Bstretching_z',function=_Bstretching_z, validators=std_validators,take_log=False)

def _normalized_LP(field,data):



print "MONKEY ON THE CAR"
