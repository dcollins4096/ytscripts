
simname = 'B20'
basedir = '/scratch1/dcollins/Paper08/%s/512'%simname
if 'frame' not in dir():
    frame = 0

in_name = "%s/RS%04d/restart%04d"%(basedir,frame,frame)
print in_name
ds = yt.load(in_name)
field_list = ['Density','x-velocity','y-velocity','z-velocity','Bx','By','Bz']
cg = ds.covering_grid(0,left_edge=[0.0,0.0,0.0],dims=[512,512,512],fields=field_list)

b02_frame = [0,10,20,30,40,50,60,70,80,87]
b2_frame = [0,10,20,30,40,50,60,63]
b20_frame = [0,10,20,30,40,50,60,69]
out_dir = "%s/RootGrids"%basedir
out_name = "%s/RootOnly_%s_512_%04d"%(out_dir,simname,frame)
out_cube = h5py.File(out_name,'w')
print out_name
for field in field_list:
    this_dataset = cg[field][:]
    this_dataset.swapaxes(0,2)
    out_cube.create_dataset(field,data=this_dataset)
out_cube.close()

