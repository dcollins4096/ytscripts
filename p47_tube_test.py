if 'ef' not in dir():
    execfile('go')
#ef('tube.py')
import tube
reload(tube)
basedir = '/Users/dcollins/scratch/Paper47_forces/'
sims = {}
sims['a02'] =  '%s/a02_dumb_test'%basedir
axis=0
coord=(0.505,0.505)
counter=0
#fields=['density','pressure','TotalEnergy','x-velocity', ]#,'x-velocity']
fields = ['Bstretching_x']
fields += ['Bx','By','Bz','magnetic_energy']
fields +=['x-velocity','y-velocity','z-velocity']
#fields += ['Ltension_x', 
#fields += ['Badvection_x']
ray_set={}
sim = 'a02'
frames=[0]
ds_list = [yt.load("%s/DD%04d/data%04d"%(sims[sim], frame,frame)) for frame in frames]
ds_list[0]
renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}
y=tube.tube(ds_list, return_ray=True,  delta=False, fields=fields, renorm=renorm,filename = "Tube_%s%s.png"%(sim,"_%s"*len(fields)%tuple(fields)), labels=fields,
            plot_args={'marker':'*'})
