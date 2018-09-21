
execfile('go')
def fake_temp(field,data):
    return data.ds.arr(np.ones_like(data['density'].v)*10,'K')

yt.add_field('temperature',function=fake_temp,units='K')
from yt.analysis_modules.ppv_cube.api import PPVCube
velocity_bounds = (-70,70,128,'code_velocity')
dims = 256
normal = [0,0,-1]
normal = 'x'
ds = yt.load('/scratch1/dcollins/Paper08/B02/256/RS0050/restart0050')
cube = PPVCube(ds, normal, "density", velocity_bounds,dims=dims, method='sum')
print "nans", np.isnan(cube[:]).sum()
