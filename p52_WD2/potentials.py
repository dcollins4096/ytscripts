

a01=taxi.taxi('p52_a01')
frame=1
class onepot():
    def __init__(self,Mi=None,Ms=None,Rb=None,Nc=None):
        self.Mi=Mi
        self.Ms=Ms
        self.Rb=Rb
        self.Nc=Nc
class potential():
    def __init__(self,taxi,verbose=False):
        self.taxi=taxi
        self.stuff={}
        self.verbose=verbose
    def read(self,frame = None):
        spherename = self.taxi.return_filename(frame)+".SphericalGravity.h5"
        if self.verbose:
            print("Opening %s"%spherename)
        if frame not in self.stuff:
            fptr=h5py.File(spherename,'r')
            try:
                Mi = fptr[ 'SphericalGravityMassInterior'][:]
                Ms = fptr['SphericalGravityMassShell'][:]
                Rb = fptr['SphericalGravityRadius'][:]
                Nc = fptr['SphericalGravityBinCount'][:]
                self.stuff[frame] = onepot(Mi=Mi,Ms=Ms,Rb=Rb,Nc=Nc)
            except:
                raise
            finally:
                fptr.close()
        return self.stuff[frame]

