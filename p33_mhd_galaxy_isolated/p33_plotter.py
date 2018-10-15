if 'ef' not in dir():
    execfile('go')

if 'framelist' not in dir():
    framelist = range(3)
if 'extra' not in dir():
    extra = -1
fieldlist = ['density','magnetic_energy']
#fieldlist += ['scaled_div_b' ]
fieldlist += ['DivB' ]
#fieldlist = ['Bx','By','Bz']
sim_base_dir = {}
sim_base_dir['e02']= 'e02_64_stock_nograckle'
sim_base_dir['e04']= 'e04_64_mhdct'
sim_base_dir['h01']= 'h01_64_stock_newCrosby'
sim_base_dir['e06']='e06_crosby_mhdct'
sim_base_dir['e07']='e07_crosby_ppm'
sim_base_dir['e08']='e08_crosby_mhdct_bsmall'
sim_base_dir['e09'] = 'e09_crosby_mhdct_bsmall_dev'
sim_base_dir['e10'] = 'e10_crosby_mhdct_bsmall_nograckle'
#bd = '/mnt/c/scratch/sciteam/dcollins/Paper33_Galaxy'
bd = '/Users/dcollins/scratch/Paper33_galaxies/'
weight_fields = {'scaled_div_b':'cell_volume'}
for sim in ['e10']: #,'e07']: #,'h01']:
    for frame in framelist:
        for ax in 'z':
            if extra > -1:
                set_name = 'ED%02d_%04d/Extra%02d_%04d'%(extra,frame,extra,frame)
                outname = 'p33_%s_ex%02d_n%04d'%(sim,extra,frame)
            else:
                set_name = 'DD%04d/DD%04d'%(frame,frame)
                outname = 'p33_%s_n%04d'%(sim,frame)
            name  = '%s/%s/%s'%(bd,sim_base_dir[sim],set_name)
            ds = yt.load(name)
            ad=ds.all_data()
            for field in fieldlist:
                stat(ad[field],field)
                stat(ad[field]-np.mean(ad[field]), field+' bar')
                proj = yt.ProjectionPlot(ds,ax,field,weight_field=weight_fields.get(field,None))
                #proj.set_zlim('density',1e-8,1e-4)
                proj.annotate_grids()
                print proj.save(outname)
