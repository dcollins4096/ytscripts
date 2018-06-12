if 'ef' not in dir():
    execfile('go')
#ef('tube.py')
import tube
reload(tube)
basedir = '/Users/dcollins/scratch/Paper49b_play/'
sims = {}
sims['a02'] =  '%s/a02_dumb_test'%basedir
sims['a03'] =  '%s/a03_transverse'%basedir
sims['r13c'] = '%s/r13c_left_fast'%basedir
sims['r14'] = '%s/r14_left_fast_test'%basedir
sims['r15'] = '%s/r15_cube_right_fast'%basedir
sims['r15b'] = '%s/r15b_test'%basedir
sims['r17'] = '%s/r17_right_python'%basedir
sims['r17b'] = '%s/r17b_right_python_test'%basedir
sims['r17c'] = '%s/r17c_check'%basedir
sims['r18'] = '%s/r18_alfven_x'%basedir
sims['r18b'] = '%s/r18b_alfven_x_fid'%basedir
sims['r19'] = '%s/r19_fast_left_rotate'%basedir
sims['r19b'] = '%s/r19b_fast_left_fid'%basedir

axis=0
coord=(0.505,0.505)
counter=0
#fields=['density','pressure','TotalEnergy','x-velocity', ]#,'x-velocity']
#fields += ['Bx','By','Bz','magnetic_energy']
#fields += ['magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_energy']
#fields +=['x-velocity','y-velocity','z-velocity']
#fields += ['Ltension_x', 
#fields += ['Badvection_x']
fields=adiabatic_mhd
#fields =['%s_p'%s for s in all_x]
ray_set={}
sim_list= ['r15','r17']
sim_list=['r15','r15b', 'r17b']
sim_list=['r17b','r17c']
sim_list=['r18b','r19']#,'r18b']
sim_list=['r19b','r19']
sim_oname = "_%s"*len(sim_list)%tuple(sim_list)
#r15=taxi.taxi('r15')
#r15b=taxi.taxi('r15b')
#r17=taxi.taxi('r17')
#flt = taxi.fleet([r15,r17,r15b])
#flt('car.load(0)')
#garr = flt('output.append( car.ds.index.grids[0])')
all_ds=[]
for this_frame in [0,50]: # [0, 34]: #,10]: # range(0,60,10):
    frames=[this_frame]

    outname = "Tubeb%s_%04d_%s.png"%(sim_oname,this_frame,"_%s"*len(fields)%tuple(fields))

    ds_list = []
    for sim in sim_list:
        ds_list += [yt.load("%s/DD%04d/data%04d"%(sims[sim], frame,frame)) for frame in frames]
    all_ds += ds_list
    print(all_ds)

    g = ds_list[0].index.grids[0]
    data=g
    renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}
    #renorm = False
    y=tube.tube(ds_list, return_ray=False,  delta=False, fields=fields, renorm=renorm,
            filename = outname,# "Tube_%s%s.png"%(sim,"_%s"*len(fields)%tuple(fields)), 
            #labels=fields,
            plot_args={'marker':'*'},legend=True, labels=sim_list)

y=tube.tube(all_ds, return_ray=False,  delta=False, fields=fields, renorm=renorm,
        filename = "TUBE_many",# "Tube_%s%s.png"%(sim,"_%s"*len(fields)%tuple(fields)), 
        #labels=fields,
        plot_args={'marker':'*'},legend=True, labels=sim_list+sim_list)
