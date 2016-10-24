if 'ef' not in dir():
    execfile('go')
#simname = 'B20'
#b02_frame = [0,10,20,30,40,50,60] #,70,80,87]
b02_frame = [70,80,87]
b2_frame = [0,10,20,30,40,50,60,63]
b20_frame = [0,10,20,30,40,50,60,69]
#b20_frame = [10,20,30,40,50,60,69]
dims = [512,512,512]
field_list = ['Density', 'PotentialField','x-velocity','y-velocity','z-velocity','Bx','By','Bz']
if 0:
    output_prefix = 'RootOnlyTake4'
    output_dir_name = 'RootGrids'
    level = 0

one_file = True
frame_list_paper08 = {'B20':b20_frame,'B2':b2_frame, 'B02':b02_frame}
if 1:
    simname = 'cubetest'
    framelist = [0]
    dims = [32,32,32]
    level = 0
    field_list = ['x-velocity']
    output_dir_name = 'RootGrids'
    output_prefix = 'OneLevel'
    one_file=True
if 0:
    print "hijack frames"
    simname = 'B20'
    b20_frame = [69]
    dims = [2048,2048,2048]
    level = 2
    print "hijack fields"
    field_list = [sys.argv[1]] #['x-velocity'] #,'y-velocity','z-velocity']
    output_dir_name = 'UnigridExtractions'
    output_prefix = 'TwoLevel'
    output_prefix += '_%s'*len(field_list)%tuple(field_list)
    one_file = False
    framelist = frame_list_paper08.get(simname,[])

if 0:
    print "hijack frames"
    b20_frame = [69]
    dims = [1024]*3
    level = 1
    print "hijack fields"
    field_list = ['x-velocity','y-velocity','z-velocity']
    output_dir_name = 'UnigridExtractions'
    output_prefix = 'OneLevel'
    one_file=True
    #output_prefix += '_%s'*len(field_list)%tuple(field_list)
    framelist = frame_list_paper08.get(simname,[])

for frame in framelist:
    basedir = '/scratch2/dcollins/Paper08/%s/512'%simname
    basedir1 = '/scratch1/dcollins/Paper08/%s/512'%simname
    basedir = '/Users/dcollins/scratch/Paper42_NewForcing/aq20_cube_test'
    basedir1 = basedir

    in_name = "%s/GravPotential/RS%04d/restart%04d"%(basedir,frame,frame)
    in_name = "%s/DD%04d/data%04d"%(basedir,frame,frame)
    out_dir = "%s/RootGrids"%basedir1
    print in_name
    ds = yt.load(in_name)

    if 0:
        print in_name
        ds = yt.load(in_name)
        #dims = [10,10,10]
        cg = ds.covering_grid(level,left_edge=[0.0,0.0,0.0],dims=dims,fields=field_list)
        out_dir = "%s/%s"%(basedir1,output_dir_name)
        if one_file:
            out_name = "%s/%s_%s_512_%04d"%(out_dir,output_prefix,simname,frame)
            out_cube = h5py.File(out_name,'w')
            print out_name
        for field in field_list:
            if not one_file:
                out_name = "%s/%s_%s_512_%04d_%s"%(out_dir,output_prefix,simname,frame,field)
                out_cube = h5py.File(out_name,'w')
                print "BUT"
                print out_name
            this_dataset = cg[field][:]
            y=this_dataset.swapaxes(0,2)
            out_cube.create_dataset(field,data=y)
            if not one_file:
                out_cube.close()
        if one_file:
            out_cube.close()

    if 1:
        grid1 = ds.index.grids[0]
        grid_direct_name = grid1.filename
        fptr_direct = h5py.File(grid_direct_name,'r')
        cube_direct = h5py.File(out_name,'r')
        s1 = fptr_direct['Grid%08d'%grid1.id]['x-velocity'][:]
        s2 = cube_direct['x-velocity'][:]
        fptr_direct.close()
        cube_direct.close()
