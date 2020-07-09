
from go import *
import taxi
import p49_QU2EB
import cmbtools
reload(taxi)
for p in ['./p42_new_driving_ak']:
    if p not in sys.path:
        sys.path += [p]
from p42_helmholtz import *
ca02 = taxi.load('ca02')
for frame in [1]:
    for axis in 'x':
        density_name= "%s/frbs/%s%04d_density_%s.fits"%(car.directory,xd,frame,axis)
        Q_name= "%s/frbs/%s%04d_Q%s.fits"%(car.directory,xd,frame,axis)
        U_name= "%s/frbs/%s%04d_U%s.fits"%(car.directory,xd,frame,axis)
        di=pyfits.open(density_name,dtype=np.double)[0].data
        qi=pyfits.open(Q_name,dtype=np.double)[0].data
        ui=pyfits.open(U_name,dtype=np.double)[0].data
        stuff=p49_QU2EB.EBfromQU(qi,ui,T=di,BoxSize=BoxSize,return_quharm=True)
