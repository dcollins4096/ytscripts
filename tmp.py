execfile('go')
def daves_pressure(field,data):
    return data['TotalEnergy']- 0.5*(data['x-velocity']**2).v
yt.add_field('daves_pressure', function= daves_pressure)
frame = 0
s1 = 'd06_ppm_zeldovich_std'
basedir = '/Users/dcollins/scratch/EnzoProjects/E12/'
sim_1 = "%s/%s"%(basedir,s1)
ds_1 = yt.load(sim_1+"/DD%04d/data%04d"%(frame,frame))
ad = ds_1.all_data()
print ad['daves_pressure'][0:100]
print ad['TotalEnergy'][0:100]
