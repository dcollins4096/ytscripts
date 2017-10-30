import re
#import h5py
import numpy as na
nar = na.array
def get_ds_name(directory, frame):
    if 0:
        name ='%s/DD%04d/DD%04d'%(directory,frame,frame)
    else:
        fname = "%s/OutputLog"%directory
        fptr = open(fname)
        lines = fptr.readlines()
        fptr.close()
        if frame >= len(lines):
            print "Dataset %d not found"%frame
            raise NotImplementedError
        setname = lines[frame].split(' ')[2]
        name = "%s/%s"%(directory,setname)
    print name
    return name
def find_grid_filename(directory, frame, grid):
    ds_name = get_ds_name(directory,frame)
    hname = "%s.hierarchy"%ds_name
    hptr = open(hname)
    grid_re = re.compile(r'^Grid = ([\d].*)')
    use_this_grid = False
    cpu_name = None
    for line in hptr:
        match = grid_re.match(line)
        if match:
            gnum= int(match.groups()[0])
            if gnum == grid:
                use_this_grid=True
                continue
        if use_this_grid:
            if line.startswith('BaryonFileName'):
                cpu_name = line.split(' ')[-1][2:-1] #skip the dots and the newline
                break
    return "%s/%s"%(directory,cpu_name)
class fake_grid():
    def __init__(self,directory, frame, grid):
        self.ds_name = get_ds_name(directory,frame)
        self.hname = "%s.hierarchy"%self.ds_name
        hptr = open(self.hname)
        grid_re = re.compile(r'^Grid = ([\d].*)')
        use_this_grid = False
        self.GridLeftEdge = None
        self.GridRightEdge = None
        self.GridLeftEdge = None
        self.GridDimension = None
        self.GridEndIndex  = None
        self.GridRank = None
        for line in hptr:
            match = grid_re.match(line)
            if match:
                gnum= int(match.groups()[0])
                if gnum == grid:
                    use_this_grid=True
                else:
                    use_this_grid=False
            if use_this_grid:
                if line.startswith('GridRank'):
                    useful_line = line[:-1].split(" ")
                    self.GridRank = int(useful_line[-1])
                    continue

                if line.startswith('GridLeftEdge'):
                    useful_line = line[:-1].split(" ")
                    self.GridLeftEdge = nar(map(float,useful_line[-(self.GridRank+1):-1]))
                    continue
                if line.startswith('GridRightEdge'):
                    useful_line = line[:-1].split(" ")
                    self.GridRightEdge = nar(map(float,useful_line[-(self.GridRank+1):-1]))
                    continue
                if line.startswith('GridStartIndex'):
                    useful_line = line[:-1].split(" ")
                    self.GridStartIndex = nar(map(int,useful_line[-(self.GridRank+1):-1]))
                    continue
                if line.startswith('GridEndIndex'):
                    useful_line = line[:-1].split(" ")
                    self.GridEndIndex = nar(map(int,useful_line[-(self.GridRank+1):-1]))
                    continue
                if line.startswith('GridDimension'):
                    useful_line = line[:-1].split(" ")
                    self.GridDimension = nar(map(int,useful_line[-(self.GridRank+1):-1]))
                    continue
        self.CellWidth = (self.GridRightEdge-self.GridLeftEdge)/(self.GridDimension-2*self.GridStartIndex)
        

def read_grid(directory, frame, grid, field):
    #pdb.set_trace()
    fname = find_grid_filename(directory,frame,grid)
    print fname
    h5ptr = h5py.File(fname)
    dset = None
    try:
        gptr = h5ptr['Grid%08d'%grid]
        dset = gptr[field][:]
    except:
        raise
    finally:
        h5ptr.close()
    return dset

def toggle(target_list,second_list):
    """swap second_list in/out of target_list"""
    for f in second_list:
        if f in target_list:
            target_list.pop( target_list.index(f) )
        else:
            target_list.append(f)
