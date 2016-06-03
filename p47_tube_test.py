if 'ef' not in dir():
    execfile('go')
#ef('tube.py')
import tube
reload(tube)
basedir = '/Users/dcollins/scratch/Paper47_forces/'
sims = {}
sims['a02'] =  '%s/a02_dumb_test'%basedir
sims['a03'] =  '%s/a03_transverse'%basedir
axis=0
coord=(0.505,0.505)
counter=0
#fields=['density','pressure','TotalEnergy','x-velocity', ]#,'x-velocity']
fields =['%s_p'%s for s in all_x]
#fields += ['Bx','By','Bz','magnetic_energy']
#fields += ['magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_energy']
#fields +=['x-velocity','y-velocity','z-velocity']
#fields += ['Ltension_x', 
#fields += ['Badvection_x']
ray_set={}
sim = 'a03'
frames=[0]
ds_list = [yt.load("%s/DD%04d/data%04d"%(sims[sim], frame,frame)) for frame in frames]
g = ds_list[0].index.grids[0]
#print g[fields[0]]
data=g
renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}
y=tube.tube(ds_list, return_ray=True,  delta=False, fields=fields, renorm=renorm,filename = "Tube_%s%s.png"%(sim,"_%s"*len(fields)%tuple(fields)), labels=fields,
            plot_args={'marker':'*'})

