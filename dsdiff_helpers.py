import re
import h5py
def get_ds_name(directory, frame):
    fname = "%s/OutputLog"%directory
    fptr = open(fname)
    lines = fptr.readlines()
    fptr.close()
    setname = lines[frame].split(' ')[2]
    return "%s/%s"%(directory,setname)
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
def find_edges(directory, frame, grid):
    ds_name = get_ds_name(directory,frame)
    hname = "%s.hierarchy"%ds_name
    hptr = open(hname)
    grid_re = re.compile(r'^Grid = ([\d].*)')
    use_this_grid = False
    GridLeftEdge = None
    GridRightEdge = None
    for line in hptr:
        match = grid_re.match(line)
        if match:
            gnum= int(match.groups()[0])
            if gnum == grid:
                use_this_grid=True
                continue
        if use_this_grid:
            if line.startswith('GridLeftEdge'):
                useful_line = line[:-1].split(" ")
                GridLeftEdge = map(float,useful_line[-4:-1])
                continue
            if line.startswith('GridRightEdge'):
                useful_line = line[:-1].split(" ")
                GridRightEdge = map(float,useful_line[-4:-1])
                continue
    return {"GridLeftEdge":GridLeftEdge,"GridRightEdge":GridRightEdge}
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
