if 'ef' not in dir():
    execfile('go')
#ef('tube.py')
import tube
reload(tube)
basedir = '/Users/dcollins/scratch/Paper33_galaxies'
basedir = '/scratch1/dcollins/Paper33_galaxies/'
basedir_35 = '/scratch1/dcollins/Paper35/'
sims = {}
sims['aj13e'] =  '%s/aj13e_square'%basedir
sims['aj13g'] =  '%s/aj13g_square_pico'%basedir
sims['aj13h'] =  '%s/aj13h_square_hancock'%basedir
sims['aj13i'] =  '%s/aj13i_square_supersonic'%basedir
sims['aj13ib'] = '%s/aj13i_repeate_supersonic_tube'%basedir_35
sims['z01'] = '%s/z01_zeldovich_fid'%basedir_35
sims['z02'] = '%s/z02_zeldovich_met'%basedir_35
sims['z03'] = '%s/z03_woodshed'%basedir_35
sims['s01'] = '%s/s01_cosmological_square'%basedir_35
sims['s03'] = '%s/s03_cosmo_square_ppm'%basedir_35
sims['s04'] = '%s/s04_ppm_large_density'%basedir_35
sims['s05'] = '%s/s05_h6_large_density'%basedir_35
sims['s06'] = '%s/s06_h0_noncosmo'%basedir_35
sims['s07'] = '%s/s07_h6_noncosmo'%basedir_35
sims['s08'] = '%s/s08_h6_hancock'%basedir_35


sim = 's01'
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
fields = ['density', 'Metal_Density','Temperature', 'TotalEnergy']# , 'metallicity']
ray_set={}
frames=range(0,10)
renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}
for sim in ['s08']:
    ds_list = []
    ds_list += [yt.load("%s/DD%04d/data%04d"%(sims[sim], frame,frame)) for frame in frames]
    fname = "Tube_%s_dve_%s.pdf"%(sim,"_%s"*len(fields)%tuple(fields))
    y=tube.tube(ds_list,   legend = True, delta=False, fields=fields,filename = fname,labels=[""]*100) #, xlim = [0.25,0.75])
#, labels=fields, plot_args={'marker':'*'})

