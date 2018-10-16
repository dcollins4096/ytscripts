import numpy as np
import mw_stuff
import uuid
import glob
import os
import pdb
import types
nar = np.array
def trans(array,dim):
    """returns the elements of array that aren't dim
    array must be a numpy array.
    trans([a,b,c],1) -> [a,c]"""
    return array[filter(lambda x: x != dim,range(len(array)) ) ]
class series():
    def __init__(self):
        self.series = []
        self.next_clump_index=0
    def next(self):
        this_clump_index = self.next_clump_index
        self.next_clump_index+=1
        if self.next_clump_index == len(self.series)+1: raise StopIteration
        return self.set[this_clump_index]
    def __iter__(self):
        self.next_clump_index=0
        return self
    def __getitem__(self,item):
        if isinstance(item,types.StringType):
            out = []
            use_clump_stuff_dict = None
            if item in ['ke/be', 'ge/all', 'ge/ke+be', 'ge/(ke+te)', 'ge/ke', 'ge/te', 'ge/be','te/be','ge/re']:
                use_clump_stuff_dict = 'Ratios'
            if item in ['ge_mw/all', 'ge_mw/ke+be', 'ge_mw/(ke+te)', 'ge_mw/ke', 'ge_mw/te', 'ge_mw/be','ge_mw/re']:
                use_clump_stuff_dict = 'Ratios_mw'
            if item in ['Gravitational', 'Kinetic', 'Rotational', 'Magnetic', 'Thermal']:
                use_clump_stuff_dict = 'Energies'
            #pdb.set_trace()
            for clump_prop in self.series:
                #Master clumps don't have 'stuff'
                try:
                    if use_clump_stuff_dict:
                        out.append(clump_prop.__dict__[use_clump_stuff_dict][item])
                    else: 
                        out.append(clump_prop.__dict__[item])
                except:
                    pass
            return nar(out)
        else:
            return self.series[item]
    def __len__(self):
        return len(self.series)
    def append(self,item):
        self.series.append(item)
class clump_stuff:
    """the clump_stuff object, sanitized to ensure code units for all quantities"""
    def __init__(self,clump,ds,Gsize_threshold=50000,clump_2d_id=-1,subclump_id=-1):
        """Attaches all quantities to the stuff object."""
        self.clump_2d_id = clump_2d_id
        self.subclump_id = subclump_id

        #For gravitational energy calculation.
        truncate = False

        if hasattr(clump,'data'):
            G = clump.data.ds['GravitationalConstant']/(4*np.pi)
            #Scalar quantites for ease
            #the mean density is almost always 1 in code units.  The extra 4 pi is for the definition of G in Enzo
            self.Time = clump.data.ds['InitialTime']
            self.TimePerFreeFall = self.Time/np.sqrt( 3*np.pi*(4*np.pi)/(32*clump.data.ds['GravitationalConstant']*1))
            data = clump.data
        elif hasattr(clump,'ds'):
            G = clump.ds['GravitationalConstant']/(4*np.pi)
            #Scalar quantites for ease
            self.Time = clump.ds['InitialTime']
            #the mean density is almost always 1 in code units.  The extra 4 pi is for the definition of G in Enzo
            data=clump
        else:
            print "nope"
            raise
        self.Time = data.ds['InitialTime']


        #set up vectors and dictionaries
        axes = np.array(['x','y','z'])
        daxes = np.array(['dx','dy','dz'])
        self.AvgColumnDensity = np.zeros(3)
        self.AvgBlos = np.zeros(3)
        truncate = False
        self.Energies = {}
        self.CenterOfMass = []

        #Check for, and move relevant zones to the right.
        self.shift = shift(data)

        #size info
        self.LeftEdge = np.array([data['x'].min(),data['y'].min(),data['z'].min()])
        self.RightEdge= np.array([data['x'].max(),data['y'].max(),data['z'].max()])
        d = self.RightEdge - self.LeftEdge
        self.MaxWidthStupid = np.sqrt(np.dot(d,d))
        self.MinWidthStupid = d.min()
        self.R = d.prod()**(1./3)*0.5
        self.R2 = 0.5*(data['cell_volume'].sum())**(1./3)

        #Projected Quantities
        for d,x in enumerate(axes):
            tx = trans(daxes,d)
            A = data[tx[0]]*data[tx[1]]
            Z = data[daxes[d]]
            self.AvgColumnDensity[d] = (A*Z*data['Density']).sum()/A.sum()
            BName = 'B%s'%('xyz'[d])
            self.AvgBlos[d] = (data['cell_mass']*data[BName]).sum()/data['cell_mass'].sum()

        #mass quantites
        self.Mass = data['cell_mass'].sum()
        for ax in ['x','y','z']:
            self.CenterOfMass.append( (data[ax]*data['cell_mass']).sum()/
                                      self.Mass)
        self.CenterOfMass = np.array(self.CenterOfMass)
        data.set_field_parameter("center",self.CenterOfMass) 
        self.mean_density = self.Mass/(data['cell_volume']).sum()


        #velocity properties 
        self.AvgVelocity = data.ds.arr([(data['cell_mass']*data[v]).sum()/data['cell_mass'].sum()
                                      for v in ['x-velocity','y-velocity','z-velocity']])
        self.VelocityUnitVector = self.AvgVelocity/np.sqrt((self.AvgVelocity**2).sum())
        data.set_field_parameter("bulk_velocity",self.AvgVelocity)
        self.VelocityDispersion =np.sqrt((data['VelocityDispersionSquared']*data['cell_mass']).sum()/self.Mass)

        #rotational properties. Require 'center' = center of mass, 'bulk velocity', both set above.
        self.AvgAngularMomentum = np.array([(data['cell_mass']*data['angular_momentum_%s'%ax]).sum()/self.Mass for ax in 'xyz'] )
        self.AvgAngularVelocity = np.array([(data['angular_momentum_%s'%ax]).sum()/self.Mass for ax in 'xyz'] )
        self.MeanOmega = np.sqrt( (self.AvgAngularVelocity**2).sum() )
        self.RotBeta = (self.MeanOmega**2*(self.R)**3)/(3*G*self.Mass)


        #test
        #self.Energies['Rotational'] = (data['cell_mass']*data['RotationalEnergy']).sum()/self.Mass
        #magnetic properties
        self.b_vol = np.zeros(3)
        self.b_mass =np.zeros(3)
        self.B_vol = np.zeros(3)
        self.B_mass= np.zeros(3)
        M = data['cell_mass']
        Mtotal = M.sum()
        V = data['cell_volume']
        Vtotal = V.sum()
        B = [data['B%s'%ax] for ax in 'xyz']
        for d in range(3):
            self.B_mass[d]= self.AvgBlos[d] #because we already did this one...
            self.B_vol[d] = (V*B[d]).sum()/Vtotal
            self.b_vol[d] =np.sqrt((V*(self.B_vol[d]-B[d].v)**2).sum()/Vtotal) 
            self.b_mass[d]=np.sqrt( (M*(self.B_mass[d]-B[d].v)**2).sum()/Mtotal)
        for f in ['cell_mass','cell_volume', 'Bx','By','Bz']:
            del data[f]

        #angles
        
        #energetic properties
        ke = (data["cell_volume"]*data["rel_kinetic_energy"]).sum()
        ke  = ke.in_units('code_velocity**2*code_mass')
        ke=ke.v
        self.Energies["Kinetic"] = ke
        be = (data["cell_volume"]*data["magnetic_energy"]).sum()
        be  = be.in_units('code_velocity**2*code_mass')
        be=be.v
        self.Energies["Magnetic"] = be
        #pdb.set_trace()
        #if data["cell_volume"].size < Gsize_threshold or Gsize_threshold < 0:
        #    ge = G*FindBindingEnergy(data["cell_mass"], data['x'],data['y'],data['z'],
        #                                                truncate, self.Energies["Kinetic"]/(G))
        #else:
        #    ge = 0
        mw_temp_name = "./MW_Grav_Temp_%s"%(uuid.uuid1())
        ge = -G*mw_stuff.run_mw_grav(mw_temp_name,data)
        
        self.Energies["Gravitational"] = ge

        #3/2nkT for monatomic gas, PV=nkt, P=cs^2 rho.
        #Assumes cs = 1
        thermal = 1.5*(data["cell_volume"]*data['Density']).sum().in_units('code_mass').v
        self.Energies["Thermal"] = thermal

        #energy ratios


        self.Ratios = {}
        if ke+be+thermal > 0.0:
            self.Ratios['ge/all']=ge/(ke+be+thermal)
        if ke+be > 0.0:
            self.Ratios['ge/ke+be']=ge/(ke+be)
        if ke > 0.0:
            self.Ratios['ge/ke']=ge/ke
        if be > 0.0:
            self.Ratios['ge/be']=ge/be
            self.Ratios['ke/be']=ke/be
            self.Ratios['te/be'] = thermal/be
        if thermal > 0.0:
            self.Ratios[ 'ge/thermal']=ge/thermal
        if ke+ge > 0.0:
            self.Ratios['ge/(ke+ge)'] = ge/(ke+ge)

        if 1: #get these later
            #if self.Energies['Rotational'] > 0.0:
            #    self.Ratios['ge/re'] = self.Energies['Gravitational']/self.Energies['Rotational']
            self.Alpha  = 5*(self.VelocityDispersion**2)/3.0*(self.R)/( self.Mass*G)
            self.Alpha2 = 5*(self.VelocityDispersion**2)/3.0*(self.R2)/( self.Mass*G)
            #self.MassToFlux = (data['MassToFluxNonCritical']*data['cell_volume']).sum()/data['cell_volume'].sum()*2*np.pi*G

        self.bound = self.Ratios['ge/all'] > 0.5
        

        #at least 2 zones wide.
        self.Valid = True
        nX = np.unique(data['x']).size
        nY = np.unique(data['y']).size
        nZ = np.unique(data['z']).size
        if (np.array([nX,nY,nZ]) <= 2).any():
            self.Valid=False
        if self.R == 0:
            self.Valid=False            

        #Need to clean up all data so i don't crash the computer.
        for k in data.keys():
            del data[k]

        mw_temp_file_list = sorted(glob.glob(mw_temp_name+"*"), reverse=True)
        for fptr in mw_temp_file_list:
            if os.path.isfile(fptr):
                os.remove(fptr)

    def __str__(self):
        out = "ge %(Gravitational).2e ke %(Kinetic).2e be %(Magnetic).2e te %(Thermal).2e"%self.Energies
        out += " ge/(all) %(ge/all).2e ge/ke %(ge/ke).2e ge/be %(ge/be).2e ke/be %(ke/be).2e"%self.Ratios
        return out

    def point_in(self,point):
        #Some error handling for old attempts.
        for i in np.arange( -self.shift[0], self.shift[0]+1):
            for j in np.arange( -self.shift[1], self.shift[1]+1):
                for k in np.arange( -self.shift[2], self.shift[2]+1):
                    delta = nar([i,j,k])* self.shift
                    if False:
                        form= "%0.2f"
                        s = "(%s, %s, %s) (%s, %s, %s) (%s, %s, %s) (%s, %s, %s)"%tuple([form]*12)
                        thrupple = self.LeftEdge.tolist()
                        thrupple+= self.RightEdge.tolist()
                        thrupple+= (point+delta).tolist()
                        thrupple+= delta.tolist()
                        print s%tuple(thrupple)
                    if (point + delta > self.LeftEdge).all() and \
                           (point + delta < self.RightEdge).all():
                        return True
        return False

def shift(clump,shiftRight = True):
    """Shifts a periodically separated clump by the domain width.
    Looks for gaps in the positions larger than max('dx'), shifts one group
    to the right (left if shiftRight=Flase) to be spatially contiguous."""
    try:
        DomainRight = clump.ds["DomainRightEdge"]
        DomainLeft =  clump.ds["DomainLeftEdge"]
        DomainWidth = DomainRight - DomainLeft
    except:
        DomainRight = clump.data.ds["DomainRightEdge"]
        DomainLeft =  clump.data.ds["DomainLeftEdge"]
        DomainWidth = DomainRight - DomainLeft
    shift = np.zeros(3)
    for i,axis in enumerate(['x','y','z']):

        dx = 'd'+axis
        nique = np.unique(clump[axis]).v
        nique.sort()
        max_dx = clump[dx].max()

        #has to be close to the edges, or 'Periodic Wrap' isn't the problem.
        if np.abs(nique.max() - DomainRight[i]) > 3*max_dx.v:
            continue
        if np.abs(nique.min() - DomainLeft[i]) > 3*max_dx.v:
            continue
        delta_x = nique[1:] - nique[0:-1]
        break_index = np.where(delta_x > max_dx.v)

        if break_index[0].size > 1:
            clump.CheckThisClump = True
        if break_index[0].size == 1:
            break_x = nique[break_index[0]]
            if shiftRight:
                all_to_shift = np.where( clump[axis] <= clump.ds.quan(break_x,'code_length') + clump[dx].min() )[0]
                clump[axis][all_to_shift] += clump.ds.arr(DomainWidth[i],'code_length')
                shift[i] = DomainWidth[i]
            else:
                all_to_shift = np.where( clump[axis] >= clump.ds.quan(break_x,'code_length') - clump[dx].min() )[0]
                clump[axis][all_to_shift] -= clump.ds.arr(DomainWidth[i],'code_length')
                shift[i] = -DomainWidth[i]
                
    try: 
        if clump.stuff:
            clump.stuff.shift = shift
    except:
        pass
    return shift
    
