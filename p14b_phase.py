
def abs_vz(field,data):
    return np.abs(data['z-velocity'])

yt.add_field('abs_vz',function=abs_vz,units='code_velocity')

if 'frame' not in dir():
    frame = 30
    simulation = 'B02'
    resolution = 256

output_prefix = 'p14b_%s_%d_n%04d'%(simulation,resolution,frame)
dirname = '/scratch1/dcollins/Paper08/%s/%d'%(simulation,resolution)
setname = "%s/RS%04d/restart%04d"%(dirname,frame,frame)
ds = yt.load(setname)

ad = ds.all_data()
phase=yt.PhasePlot(ad,'density','abs_vz','cell_mass',weight_field=None)
print phase.save(output_prefix)


