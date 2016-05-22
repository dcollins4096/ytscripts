
import types
import numpy as na
import pdb
nar = na.array


class subset():
    def __init__(self,clumpset):
        self.clumpset = clumpset
        self.name='empty'
        self.set = []

        #Two kludges
        self.all = self
        self.most_single = self
        self.next_clump_index=0
    def next(self):
        this_clump_index = self.next_clump_index
        self.next_clump_index+=1
        if self.next_clump_index == len(self.set)+1: raise StopIteration
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
            for clump_prop in self.clumpset:
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
            return self.set[item]
    def __len__(self):
        return len(self.set)
class all(subset):
    def __init__(self,clumpset,skip_first=True):
        """Generate a flat list of all clumps.  
        For clumps sharing a large parent, skip_first ignores the first clump
        which is often the entire domain."""
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'all'
        first_index = 0
        if skip_first:
            first_index = 1

        for tree in self.clumpset.all_for_master:
            self.set.extend(tree[first_index:])
        self.set=nar(self.set)

class all_valid(subset):
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.name='all_valid'
        self.set = clumpset.all[ clumpset.all['Valid'] ]

class most_single(subset):            
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'most_single'
        for c in clumpset.all:
            if c.parent == None:
                continue
            if len(c.parent.children) == 1:
                continue
            if c.children is not None:
                if len(c.children) > 1:
                    continue
            self.set.append(c)
        self.set = nar(self.set)
class energetically_bound(subset):            
    """All single clumps where (ke+te+be)/ge < 2"""
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'energetically_bound'

        energetic_list=na.where( nar([1./(x['ge/all']) for x in clumpset.most_single['Ratios']]) < 2)[0]
        self.set = clumpset.most_single[energetic_list]
class single_unbound(subset):
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'single_unbound'
        energetic_list=na.where( nar([1./(x['ge/all']) for x in clumpset.most_single['Ratios']]) > 2)[0]
        self.set = clumpset.most_single[energetic_list]
class energetically_bound_2(subset):
    """All clumps where (ke+te+be)/ge < 2"""
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'energetically_bound_2'

        energetic_list=na.where( nar([1./(x['ge/all']) for x in clumpset.all['Ratios']]) < 2)[0]
        self.set = clumpset.all[energetic_list]
def has_bound_children(clump):
    """If any """
    for c in clump.children:
        if c.stuff.Ratios['ge/all'] > 0.5:
            return True
        if has_bound_children(c):
            return True
    return False
class energetically_bound_3(subset):
    """All clumps where (ke+te+be)/ge < 2 that don't have bound children."""
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'energetically_bound_3'

        energetic_list=na.where( nar([1./(x['ge/all']) for x in clumpset.all['Ratios']]) < 2)[0]
        set_all = clumpset.all[energetic_list]
        if len(set_all) == 0:
            self.set = set_all
            return
        #cull this set for children that are bound
        keep = nar([True]* len(set_all))
        for i,c in enumerate(set_all):
            if has_bound_children(c):
                keep[i] = False
        self.set = set_all[keep]
class alpha_bound(subset):            
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'alpha_bound'

        alpha_list=na.where( clumpset.most_single['Alpha'] < 2)[0]
        self.set = clumpset.most_single[alpha_list]
class violated(subset):            
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'violated'

        violated_list=na.where( clumpset.most_single['ViolatedMass'] != 0 )[0]
        self.set = clumpset.most_single[violated_list]

class globally_violated(subset):            
    """All clumps with non-zero ViolatedMass and refined to the max level."""
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'globally_violated'

        is_violated =clumpset.most_single['ViolatedMass'] != 0 
        max_level = clumpset.most_single[0].data.pf.h.max_level
        is_max_level = [c['GridLevel'].max() == max_level for c in clumpset.most_single]
        violated_list = na.logical_and( is_max_level, is_violated)
        self.set = clumpset.most_single[violated_list]
class not_so_violated_bound(subset):
    """Energetically bound, with M_{v}/M < 0.3"""
    def __init__(self,clumpset):
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'not_so_violated_bound'

        energetic_list= nar([1./(x['ge/all']) for x in clumpset.most_single['Ratios']]) < 2
        violated_list = clumpset.most_single['ViolatedMass']/clumpset.most_single['MassGram'] < 0.3
        both = na.logical_and( energetic_list,violated_list)
        self.set = clumpset.most_single[both]
#specifically for paper 5
class small_unbound(subset):
    def __init__(self,clumpset):
        """ unbound, r<4e-3"""
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'small_unbound'
        energetic_list=nar([1./(x['ge/all']) for x in clumpset.most_single['Ratios']]) > 2
        small = clumpset.most_single['R'] < 4e-3
        self.set = clumpset.most_single[ na.logical_and( energetic_list, small)]
class large_unbound(subset):
    def __init__(self,clumpset):
        """ unbound, r>4e-3"""
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'large_unbound'
        energetic_list=nar([1./(x['ge/all']) for x in clumpset.most_single['Ratios']]) > 2
        small = clumpset.most_single['R'] > 4e-3
        self.set = clumpset.most_single[ na.logical_and( energetic_list, small)]

class biggest(subset):
    def __init__(self,clumpset):
        """ Lowest level material: master list and the first children."""
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'biggest'
        for c in clumpset.master_list:
            self.set.append(c)
            for child in c.children:
                self.set.append(child)

class TooBigToMeasure(subset):
    def __init__(self,clumpset):
        """Where Gravitational is zero"""
        subset.__init__(self,clumpset)
        self.set = []
        self.name = 'too_big'
        self.set = clumpset.all[ (clumpset.all['Gravitational'] == 0 ) ]
