#from yt.mods import *
import yt
from math import *
import numpy as np
nar = np.array
import pdb
#from collections import defaultdict
#import uber
#from davetools import *
#My version of Clump is used in order to preserve the 'stuff' pickling/unpickling.
#import DaveClump 
#import clump_subset
#from yt.utilities import physical_constants as units
from yt.utilities.data_point_utilities import FindBindingEnergy

"""2010/12/06 Removed the following:
    tabs(level)  :can be replaced with simpler code.
    multi_cores  :not presently in use.  Maybe more useful later?
    binary_iterations from all clump_set stuff.
    """

#for recursive functions.
counter=0


def clump_list_sort(clump_list):
    minDensity = [c['Density'].min() for c in clump_list]
    sortedargs = np.argsort(minDensity)
    list = np.array(clump_list[:])[sortedargs]
    reverse = range(list.size-1,-1,-1)
    return list[reverse]
                    

def find_clumps(clump, min, max, d_clump, method='mult'):
    print "Finding clumps: min: %e, max: %e, step: %f" % (min, max, d_clump)
    if min >= max: return
    clump.find_children(min)

    if method == 'mult':
        next_level = min*d_clump
    elif method == 'add':
        next_level = min + d_clump
    print "Finding clumps: next level min %e next_level %e" % (min, next_level)
    if (len(clump.children) == 1):
        find_clumps(clump, next_level, max, d_clump,method)
    elif (len(clump.children) > 0):

        these_children = []
        print "Investigating %d children." % len(clump.children)
        for child in clump.children:
            find_clumps(child, min*d_clump, max, d_clump,method)
            if ((child.children is not None) and (len(child.children) > 0)):
                these_children.append(child)
            elif (child._isValid()):
                these_children.append(child)
            else:
                print "Eliminating invalid, childless clump with %d cells." % len(child.data["CellMassMsun"])
        if (len(these_children) > 1):
            print "%d of %d children survived." % (len(these_children),len(clump.children))
            clump.children = these_children
        elif (len(these_children) == 1):
            print "%d of %d children survived, linking its children to parent." % (len(these_children),len(clump.children))
            clump.children = these_children[0].children
        else:
            print "%d of %d children survived, erasing children." % (len(these_children),len(clump.children))
            clump.children = []

            
def return_bottom_clumps(clump,dbg=0):
    global counter
    counter = 0
    list = []
    level = 0
    recursive_bottom_clumps(clump,list,dbg,level)
    return list
def recursive_bottom_clumps(clump,clump_list, dbg = 0,level=0):
    """Loops over a list of clumps (clumps) and fills clump_list with the bottom most.
    Recursive. Prints the level and the number of cores to the screen."""

    global counter
    if dbg > 0:
        print "\t"*level, "l =",level, "n_core",counter

    if ((clump.children is None) or (len(clump.children) == 0)):
        counter += 1
        clump_list.append( clump )
    else:
        for child in clump.children:
            recursive_bottom_clumps(child,clump_list,dbg=dbg,level=level+1)

def recursive_all_clumps(clump,list,level,parentnumber,numberClumps):
    global counter
    counter += 1
    if numberClumps:
        clump.number = counter
        clump.parentnumber = parentnumber
        clump.level = level
    counter += 1
    list.append(clump)
    if clump.children != None:
        for child in clump.children:
            x = recursive_all_clumps(child,list,level+1,clump.number,numberClumps)
    return list

def return_all_clumps(clump,numberClumps = True):
    global counter
    counter = 0
    list = []
    level = 0
    clump.level = level
    parentnumber=-1
    recursive_all_clumps(clump,list,level,parentnumber,numberClumps)
    return list


import fPickle
import glob

class clump_set:
    def __init__(self,data_object=None,field='Density',
                 step_size = None, c_min = None, c_max = None,n_steps = None,
                 function = 'True',input_master=None,non_recursive=False,huge_master=True,method='mult'):
        """
        This is a wrapper for the clump finder in Clump.py.
        Arguments:
        data_object: object over which contouring is performed (region or sphere).
        field: data field over which contours are made (example: "Density" or "AveragedDensity").
        step: contouring stepsize.  The field minimum is multiplied by this value each round of 
        the clump finding.
        *input_master* overrides the clump finder and generates a clump_set from the input master_list
        *huge_master* skips the first clump when making subsets, since this is frequently the entire domain 
            and not of interest.
        """

        if input_master == None:
            self.bound = None
            
            print n_steps
            if method == 'mult':
                if c_min == None:
                    c_min = 10**floor(log10(data_object[field].min()))
                if c_max == None:
                    c_max = 10**floor(log10(data_object[field].max())+1)
                if n_steps != None:
                    step_size = 10**( (np.log10( c_max ) - np.log10( c_min ) )/n_steps )
                    #step_size = ( c_max - c_min )/n_steps 
                else:
                    n_steps = (np.log10(c_max) - np.log10(c_min))/np.log10(step_size)
                if step_size == None:
                    print "please specifiy either n_steps or step_size."
                    return None
            elif method == 'add':
                if c_min == None:
                    c_min = floor(data_object[field].min())
                if c_max == None:
                    c_max = floor(data_object[field].max())+1
                if n_steps != None:
                    step_size = (c_max - c_min ) /n_steps
                    #step_size = ( c_max - c_min )/n_steps 
                else:
                    n_steps = (c_max - c_min ) / step_size
                if step_size == None:
                    print "please specifiy either n_steps or step_size."
                    return None

            
            self.c_min = c_min
            self.c_max = c_max
            self.step_size = step_size
            self.n_steps = n_steps

            self.levels = c_min*step_size**np.arange(n_steps)

            self.function = "np.logical_and(self['Density'].size > 27, "
            self.function += "(np.array([np.unique(self['x']).size,np.unique(self['y']).size,np.unique(self['z']).size])"
            self.function += " >= 3).all())"

            master_list = []
            counter=0
            for data in data_object:
                counter+=1
                print "Finding clumps, %d/%d"%(counter,len(data_object))
                master_list.append( DaveClump.Clump(data, None, field,
                                                cached_fields = defaultdict(lambda: dict()),
                                                function = self.function))
            
                if method == 'non_recursive':
                    master_list[-1].find_children(c_min,c_max)
                else:
                    find_clumps(master_list[-1], c_min, c_max, step_size,method)

            #find_dendrogram(master_clump, c_min, c_max, step_size)


            self.master_list = master_list
        else:
            self.master_list = input_master
        
            
        #self.finest = np.array(return_bottom_clumps(self.master))
        self.all_for_master = []
        for master in self.master_list:
            
            self.all_for_master.append(return_all_clumps(master))

        self.huge_master = huge_master
        self.units = None
        #the need to skip, or not, the first subset means that the "all"
        #subset must be treated alone.  All other subsets may use the add_subset mechanism.
        #self.add_subset(clump_subset.all)
        self.all = clump_subset.all(self,skip_first=huge_master)
        self.add_subset(clump_subset.most_single)
        self.number_tracks()
        
    def set_units(self,density=None,velocity=None,length=None,temperature=None):
        """sets multipliers for derived quantities. We'll need to come back to this. """

        pass

    def add_subset(self,subset):
        subset_instance=subset(self)
        if self.__dict__.has_key(subset_instance.name):
            print "This clump list has something by the name "+subset_instance.name
        else:
            self.__dict__[subset_instance.name]=subset_instance
    def flagged(self,flag,quantity=None):
        if quantity is not None:
            all_q = self.all[quantity]
        output = []
        for n, cl in enumerate(self.all):

            if hasattr(cl,'sets'):
                if flag in cl.sets:
                    if quantity:
                        output.append( all_q[n-1])
                    else:
                        output.append(cl)
        return nar(output)


    def number_tracks( self ):
        for n, c in enumerate(self.all):
            c.track = n
    def return_tracks(self):
        return np.array( [c.track for c in self.all] )

    def select_tracks(self, track_list):
        out = []
        track_listing = self.return_tracks()
        for track in track_list:
            out += list(self.all[ track_listing == track ]) 
        return out

    def stuff(self, size_threshold=50000,clobber=False):
        """Compute clump_stuff for each clump in the hierarchy.
        size_threshold is the maximum size for which binding energy will be computed.
        Also create self.valid array of clumps for which clump_stuff.valid is true.
        *clobber* will eliminate existing stuff instances"""
        self.valid = []
        counter=0
        for a in self.all:
            counter += 1
            print "stuff on %d/%d"%(counter,len(self.all))
            if not hasattr(a,'stuff') or clobber == True:
                a.stuff = clump_stuff_code(a,size_threshold)
            else: 
                print "    extant.  To remake, clobber = True"
            try:
                if a.stuff.Valid:
                    self.valid.append(a)
            except:
                pass
        self.valid = np.array(self.valid)
    def save(self,filename):
        """Saves just the master list of clump_set to file called *filename*.
        Very simple function, ensures a match with clump_stuff.load"""
        fPickle.bdump(self.master_list,filename)

def load(filename):
    """Load master list from *filename*, populate clump_set"""
    full_list = fPickle.bload(filename)
    master_list = []
    """I'm not sure if this is right."""
    huge_master = False 
    try:
        len(full_list[0])
        pf = full_list[0][0]
        for master in full_list:
            master_list.append(master[1])
            set = clump_set(input_master = master_list,huge_master=huge_master)
        set.pf = pf            
    except:
        try:
            pf = full_list[0]
            master_list = [full_list[1]]
            set = clump_set(input_master = master_list,huge_master=huge_master)
            set.pf = pf
        except:
            set = full_list

    #the existance of this pf is maybe unreliable right now,
    #but it definitely needs to be saved.

    return set

class clump_stuff_code:
    """the clump_stuff object, sanitized to ensure code units for all quantities"""
    def __init__(self,clump,ds,Gsize_threshold=50000):
        """Attaches all quantities to the stuff object."""

        #For gravitational energy calculation.
        truncate = False

        self.Time = clump.data.pf['InitialTime']
        self.TimePerFreeFall = self.Time/np.sqrt( 3*np.pi/(32*G))
        if hasattr(clump,'data'):
            G = clump.data.pf['GravitationalConstant']/(4*np.pi)
            #Scalar quantites for ease
            #the mean density is almost always 1 in code units.  The extra 4 pi is for the definition of G in Enzo
            self.TimePerFreeFall = self.Time/np.sqrt( 3*np.pi*(4*np.pi)/(32*clump.data.pf['GravitationalConstant']*1))
            convert_to_cm = clump.data.convert('cm')
            data = clump.data
        elif hasattr(clump,'pf'):
            G = clump.pf['GravitationalConstant']/(4*np.pi)
            #Scalar quantites for ease
            self.Time = clump.pf['InitialTime']
            #the mean density is almost always 1 in code units.  The extra 4 pi is for the definition of G in Enzo
            convert_to_cm = clump.convert('cm')
            data=clump
        else:
            print "nope"
            raise


        #set up vectors and dictionaries
        axes = np.array(['x','y','z'])
        daxes = np.array(['dx','dy','dz'])
        self.AvgColumnDensity = np.zeros(3)
        self.AvgBlos = np.zeros(3)
        truncate = False
        self.Energies = {}
        self.CenterOfMass = []

        #Check for, and move relevant zones to the right.
        self.shift = shift(clump)

        #size info
        self.LeftEdge = np.array([clump['x'].min(),clump['y'].min(),clump['z'].min()])
        self.RightEdge= np.array([clump['x'].max(),clump['y'].max(),clump['z'].max()])
        d = self.RightEdge - self.LeftEdge
        self.MaxWidthStupid = np.sqrt(np.dot(d,d))
        self.MinWidthStupid = d.min()
        self.R = d.prod()**(1./3)*0.5
        self.R2 = 0.5*(clump['CellVolume'].sum())**(1./3)

        #Projected Quantities
        for d,x in enumerate(axes):
            tx = trans(daxes,d)
            A = clump[tx[0]]*clump[tx[1]]
            Z = clump[daxes[d]]*convert_to_cm
            self.AvgColumnDensity[d] = (A*Z*clump['Density']).sum()/A.sum()
            self.AvgBlos[d] = (clump['CellMass']*clump['MagneticField_C_%1d'%(d+1)]).sum()/clump['CellMass'].sum()

        #mass quantites
        self.Mass = clump['CellMass'].sum()
        for ax in ['x','y','z']:
            self.CenterOfMass.append( (clump[ax]*clump['CellMass']).sum()/
                                      self.Mass)
        self.CenterOfMass = np.array(self.CenterOfMass)
        data.set_field_parameter("center",self.CenterOfMass) 
        self.mean_density = self.Mass/(clump['CellVolume']).sum()


        #velocity properties 
        self.AvgVelocity = np.array([(clump['CellMass']*clump[v]).sum()/clump['CellMass'].sum()
                                      for v in ['x-velocity','y-velocity','z-velocity']])
        data.set_field_parameter("bulk_velocity",self.AvgVelocity)
        self.VelocityNorm = self.AvgVelocity/np.sqrt((self.AvgVelocity**2).sum())
        self.VelocityDispersion =data.quantities['rmsVelocityDispersionMass']()

        #rotational properties. Require 'center' = center of mass, 'bulk velocity', both set above.
        self.AvgAngularMomentum = np.array([(clump['CellMass']*clump['AngularMomentum'][i,:]).sum()/self.Mass for i in range(3)] )
        self.AvgAngularVelocity = np.array([(clump['CellMass']*clump['AngularVelocity'][i,:]).sum()/self.Mass for i in range(3)] )
        self.MeanOmega = np.sqrt( (self.AvgAngularVelocity**2).sum() )
        self.RotBeta = (self.MeanOmega**2*(self.R)**3)/(3*G*self.Mass)


        #test
        self.Energies['Rotational'] = (clump['CellMass']*clump['RotationalEnergy']).sum()/self.Mass
        #magnetic properties
        self.b_vol = np.zeros(3)
        self.b_mass =np.zeros(3)
        self.B_vol = np.zeros(3)
        self.B_mass= np.zeros(3)
        M = clump['CellMass']
        Mtotal = M.sum()
        V = clump['CellVolume']
        Vtotal = V.sum()
        B = [clump['MagneticField_C_%1d'%(d)] for d in range(1,4)]
        for d in range(3):
            self.B_mass[d]= self.AvgBlos[d] #because we already did this one...
            self.B_vol[d] = (V*B[d]).sum()/Vtotal
            self.b_vol[d] =np.sqrt((V*(self.B_vol[d]-B[d])**2).sum()/Vtotal) 
            self.b_mass[d]=np.sqrt( (M*(self.B_mass[d]-B[d])**2).sum()/Mtotal)
        for f in ['CellMass','CellVolume']+['MagneticField_C_%1d'%(d) for d in range(1,4)]: 
            del data[f]

        #angles
        
        #energetic properties
        ke = (clump["CellVolume"]*clump["RelKineticEnergy"]).sum()
        self.Energies["Kinetic"] = ke
        be = (clump["CellVolume"]*clump["MagneticEnergy"]).sum()
        self.Energies["Magnetic"] = be
        #pdb.set_trace()
        if clump["CellVolume"].size < size_threshold or size_threshold < 0:
            ge = G*FindBindingEnergy(clump["CellMass"], clump['x'],clump['y'],clump['z'],
                                                        truncate, self.Energies["Kinetic"]/(G))
        else:
            ge = 0
        self.Energies["Gravitational"] = ge

        #3/2nkT for monatomic gas, PV=nkt, P=cs^2 rho.
        #Assumes cs = 1
        thermal = 1.5*(clump["CellVolume"]*clump['Density']).sum()
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
        if self.Energies['Rotational'] > 0.0:
            self.Ratios['ge/re'] = self.Energies['Gravitational']/self.Energies['Rotational']

        self.Alpha  = 5*(self.VelocityDispersion**2)/3.0*(self.R)/( self.Mass*G)
        self.Alpha2 = 5*(self.VelocityDispersion**2)/3.0*(self.R2)/( self.Mass*G)
        self.MassToFlux = (clump['MassToFluxNonCritical']*clump['CellVolume']).sum()/clump['CellVolume'].sum()*2*np.pi*G

        self.bound = self.Ratios['ge/all'] > 0.5
        

        #at least 2 zones wide.
        self.Valid = True
        nX = np.unique(clump['x']).size
        nY = np.unique(clump['y']).size
        nZ = np.unique(clump['z']).size
        if (np.array([nX,nY,nZ]) <= 2).any():
            self.Valid=False
        if self.R == 0:
            self.Valid=False            

        #Need to clean up all data so i don't crash the computer.
        for k in data.keys():
            del data[k]

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
        DomainRight = clump.pf["DomainRightEdge"]
        DomainLeft =  clump.pf["DomainLeftEdge"]
        DomainWidth = DomainRight - DomainLeft
    except:
        DomainRight = clump.data.pf["DomainRightEdge"]
        DomainLeft =  clump.data.pf["DomainLeftEdge"]
        DomainWidth = DomainRight - DomainLeft
    shift = np.zeros(3)
    for i,axis in enumerate(['x','y','z']):

        dx = 'd'+axis
        nique = np.unique(clump[axis])
        nique.sort()
        max_dx = clump[dx].max()

        #has to be close to the edges, or 'Periodic Wrap' isn't the problem.
        if np.abs(nique.max() - DomainRight[i]) > 3*max_dx:
            continue
        if np.abs(nique.min() - DomainLeft[i]) > 3*max_dx:
            continue
        delta_x = nique[1:] - nique[0:-1]
        break_index = np.where(delta_x > max_dx )

        if break_index[0].size > 1:
            clump.CheckThisClump = True
        if break_index[0].size == 1:
            break_x = nique[break_index[0]]
            if shiftRight:
                all_to_shift = np.where( clump[axis] <= break_x + clump[dx].min() )[0]
                clump[axis][all_to_shift] += DomainWidth[i]
                shift[i] = DomainWidth[i]
            else:
                all_to_shift = np.where( clump[axis] >= break_x - clump[dx].min() )[0]
                clump[axis][all_to_shift] -= DomainWidth[i]
                shift[i] = -DomainWidth[i]
                
    try: 
        if clump.stuff:
            clump.stuff.shift = shift
    except:
        pass
    return shift
    
