

#g = 
fft_func=np.fft.fftn
aq24 = taxi.taxi(dir='/Users/dcollins/scratch/Paper42_NewForcing/aq24_sto_b',name='p42_aq24')
aq24.fill(20)
g=aq24.ds.index.grids[0]
v1 = g['x-velocity']
v2 = g['y-velocity']
v3 = g['z-velocity']
v1hat = fft_func(v1)
v2hat = fft_func(v2)
v3hat = fft_func(v3)
if 0:
    print ((v1*v1)).sum()/(v1hat*np.conj(v1hat)).sum()*v1.size
    v1hathat = np.fft.ifftn(v1hat)

div_fac = 2.0
dx = 1./32
ftype = 'gas'
sl_center = slice(1, -1, None)
sl_left = slice(None, -2, None)
sl_right = slice(2, None, None)
div_fac = 2.0
ftype = 'gas'

def d_dx(arr):
    ds = 2*1./32 #div_fac * just_one(data["index", "dx"])
    f  = arr[sl_right,1:-1,1:-1]/ds
    f -= arr[sl_left ,1:-1,1:-1]/ds
    new_field = np.zeros_like(arr)
    new_field[1:-1,1:-1,1:-1] = f
    return new_field
def d_dy(arr):
    ds = 2*1./32 #div_fac * just_one(data["index", "dx"])
    f  = arr[1:-1,sl_right,1:-1]/ds
    f -= arr[1:-1,sl_left ,1:-1]/ds
    new_field = np.zeros_like(arr)
    new_field[1:-1,1:-1,1:-1] = f
    return new_field

def d_dz(arr):
    ds = 2*1./32 #div_fac * just_one(data["index", "dx"])
    f  = arr[1:-1,1:-1,sl_right]/ds
    f -= arr[1:-1,1:-1,sl_left ]/ds
    new_field = np.zeros_like(arr)
    new_field[1:-1,1:-1,1:-1] = f
    return new_field

divv=d_dx(v1)+ d_dy(v2)+ d_dz(v3)

nx = 1.shape
kvec = np.ogrid[0:nx[0],0:nx[1],0:nx[2]]
NormK = kvec[0]**2+kvec[1]**2+kvec[2]**2
NormK[0,0,0]=1.0


