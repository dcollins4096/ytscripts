import yt
import numpy as np
def specific_ke(field,data):
    vx = data['x-velocity'].in_units('code_velocity').v
    vy = data['y-velocity'].in_units('code_velocity').v
    vz = data['z-velocity'].in_units('code_velocity').v
    return 0.5*(vx*vx+vy*vy+vz*vz)
yt.add_field('specific_ke',specific_ke,take_log=True)
def ge_dve(field,data):
    try:
        return data[('enzo','GasEnergy')].v
    except:
        vx = data['x-velocity'].in_units('code_velocity').v
        vy = data['y-velocity'].in_units('code_velocity').v
        vz = data['z-velocity'].in_units('code_velocity').v
        return data[('enzo','TotalEnergy')].v-0.5*(vx*vx+vy*vy+vz*vz)

def p_dave(field,data):
    vx = data['x-velocity'].in_units('code_velocity').v
    vy = data['y-velocity'].in_units('code_velocity').v
    vz = data['z-velocity'].in_units('code_velocity').v
    p =data[('enzo','TotalEnergy')].v-0.5*(vx*vx+vy*vy+vz*vz)
    p*= data[('enzo','Density')]/(data.ds['Gamma']-1)
    return p
yt.add_field('p_dave',p_dave)

def s_dave(field,data):
    p = data['p_dave']
    d = data['density'].in_units('code_density').v
    gamma = data.ds['Gamma']
    p_over_rho = p/d**gamma
    if 0:
        try:
            min_p = (p_over_rho[ p_over_rho > 0]).min()
            p_over_rho[ p_over_rho <= 0] = min_p
        except:
            pass
    s = p_over_rho # -np.log(p_over_rho)
    return s
yt.add_field('s_dave',s_dave,take_log=False)

yt.add_field('ge',ge_dve, take_log=False)
def nolog_ke(field,data):
    return data['kinetic_energy'].v
yt.add_field('ke_nolog', nolog_ke, take_log=False)
def temp_vs_ge(field,data):
    return data[('enzo','Temperature')].v/data['ge']
yt.add_field('tvg',temp_vs_ge)
def te_enzo(field,data):
    return data[('enzo','TotalEnergy')].v
yt.add_field('te_enzo',te_enzo,take_log=False)
def te_enzo_b(field,data):
    return data[('enzo','TotalEnergy')].v
yt.add_field('te_enzo_log',te_enzo,take_log=True)

class get_the_stuff():
    def __init__(self,fname_in):
        fname = "%s.cpu0000"%(fname_in)
        print fname
        fptr = h5py.File(fname)
        lots=True
        mag=False
        try:
#        if True:
            grid = fptr['Grid%08d'%(1)]
            self.temperature = grid['Temperature'][:]
            ge_disk = None
            if lots:
                self.density = grid['Density'][:]
                self.vx = grid['x-velocity'][:]
                self.vy = grid['y-velocity'][:]
                self.vz = grid['z-velocity'][:]
                self.te = grid['TotalEnergy'][:]
                if 'GasEnergy' in grid.keys(): 
                    print "I have gas"
                    self.ge_disk = grid['GasEnergy'][:]
                    stat(self.ge_disk,"Gas Energy, Disk")
                else:
                    print "Gas-X.  "
                self.ke = 0.5*(self.vx*self.vx+self.vy*self.vy+self.vz*self.vz)
                self.ge_find = self.te - self.ke
                stat(self.ge_find,"Gas energy, Derived.")
            if mag:
                self.bx = grid['Bx'][:]
                self.by = grid['By'][:]
                self.bz = grid['Bz'][:]
                self.be = 0.5*(self.bx*self.bx+self.by*self.by+self.bz*self.bz)
                stat( be, "BE")

        except:
            raise
        finally:
            fptr.close()
