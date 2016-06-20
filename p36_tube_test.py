if 'ef' not in dir():
    execfile('go')
#ef('tube.py')
import tube
reload(tube)
sim_4 = '/Users/dcollins/scratch/Paper36/t04_pressure_hydro6'; simname='t04'
sim_3 = '/Users/dcollins/scratch/Paper36/t03_pressure'; simname='t03'
sim_5 = '/Users/dcollins/scratch/Paper36/t05_sine_fast'; simname='t05'
#ds_3 = yt.load('%s/DD%04d/data%04d'%(sim_0,n0,n0))
if 0:
    n0 = 0
    n1 = 10
    if 'ds0' not in dir() or True:
        ds0 = yt.load('%s/DD%04d/data%04d'%(sim,n0,n0))
        ds1 = yt.load('%s/DD%04d/data%04d'%(sim,n1,n1))
axis=0
coord=(0.505,0.505)
counter=0
#fields=['density','pressure','TotalEnergy','x-velocity', ]#,'x-velocity']
#fields += ['By','Bz', 'velocity_y','velocity_z']
ray_set={}
#ds_list=[ds0,ds1]
subset = [0,2,4]+[10] + [100]
ds_list_mhd = [yt.load('%s/DD%04d/mhd%04d'%(sim_4,n1,n1)) for n1 in subset]
ad_mhd0 = ds_list_mhd[0].all_data()
ds_list_hydro = [yt.load('%s/DD%04d/data%04d'%(sim_3,n1,n1)) for n1 in subset]
ds_list_sine = [yt.load('%s/DD%04d/data%04d'%(sim_5,n1,n1)) for n1 in subset]
ds_list = ds_list_sine #[ds_list_mhd[0], ds_list_hydro[0]]
renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}
if 0:
    ds0 = ds_list[0]
    ds1 = ds_list[1]
    for ds in ds_list:
        ray_set[ds]= ds.ortho_ray(axis,coord) 
    ray_mhd=ray_set[ds0]
    ray_hydro=ray_set[ds1]
    simname = "both_ics"
simname = 't05_sine'
#x=tube.tube(ds_list,      delta=True,  fields=fields, renorm=renorm,filename = simname+".png",legend=True)
#y=tube.tube(ds_list_mhd,  delta=true, fields=fields, renorm=renorm,filename = "mhd"+".png")
y=tube.tube(ds_list_sine,  delta=True, fields=fields, renorm=renorm,filename = "sine"+".png")
#y=tube.tube(ds_list_hydro,delta=True,  fields=fields, renorm=renorm,filename = "hydro"+".png")
if 0:
    for i,field in enumerate(fields):
        for n_ds,ds in enumerate(ds_list):
            counter += 1
            this_ray = ray_set[ds]
            this_x = this_ray['xyz'[axis]]
            this_y = this_ray[field]
            stat(this_y)

#needs: ensure list
