if 'ef' not in dir():
    execfile('go')
#ef('tube.py')
import tube
reload(tube)
basedir = '/Users/dcollins/scratch/Paper33_galaxies'
basedir = '/scratch1/dcollins/Paper33_galaxies/'
sims = {}
sims['aj13e'] =  '%s/aj13e_square'%basedir
sims['aj13g'] =  '%s/aj13g_square_pico'%basedir
sims['aj13h'] =  '%s/aj13h_square_hancock'%basedir
sims['aj13i'] =  '%s/aj13i_square_supersonic'%basedir
sim = 'aj13i'
axis=0
coord=(0.505,0.505)
counter=0
#fields=['density','pressure','TotalEnergy','x-velocity', ]#,'x-velocity']
#fields =['%s_p'%s for s in all_x]
#fields += ['Bx','By','Bz','magnetic_energy']
#fields += ['magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_energy']
#fields +=['x-velocity','y-velocity','z-velocity']
#fields += ['Ltension_x', 
#fields += ['Badvection_x']
fields = ['density','Metal_Density', 'metallicity']
ray_set={}
frames=[0,3,10] #range(10) #[0,1,5,10] #,2,5] #range(5) #[0,10] #range(30)
renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}

ds_list = [yt.load("%s/DD%04d/data%04d"%(sims[sim], frame,frame)) for frame in frames]
y=tube.tube(ds_list,   delta=False, fields=fields,filename = "Tube_%s%s.png"%(sim,"_%s"*len(fields)%tuple(fields)),labels=[""]*100)
#, labels=fields, plot_args={'marker':'*'})

