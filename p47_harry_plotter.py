if 'ef' not in dir():
    execfile('go')

import dsdiff_helpers
import dsdiff
from mpl_toolkits.axes_grid1 import AxesGrid
reload(dsdiff)
basedir = '/Users/dcollins/scratch/Paper47_forces/'
simdict={}
simdict['a01'] = 'a01'
simdict['a04'] = 'a04_collision_weak_field'
simdict['a05'] = 'a05_collision_no_field'
simdict['a06'] = 'a06_collision_strong_field'
simdict['a07'] = 'a07_strong_field_medium_density'
simdict['a07b'] = 'a07b_strong_periodic'
simdict['a08'] = 'a08_oblique'
simlist = ['a04', 'a07']
if 'extra' not in dir():
    extra = -1
if 'framelist' not in dir():
    framelist = [0, 100] #range(90,101)
    framelist = [0]

fieldlist = []
#fieldlist = ['density']
#fieldlist += ['Bx']
#fieldlist += ['magnetic_energy']
#fieldlist += ['By','Bz']
#fieldlist = [('enzo','DivB')]
#fieldlist += ['density']
#fieldlist += [('enzo','TotalEnergy')]
if 0:
    fieldlist = ['Badvection_x']
    fieldlist += ['Badvection_y']
    fieldlist += ['Badvection_z']
if 0:
    fieldlist += ['%s-velocity'%s for s in 'xyz']
if 0:
    fieldlist = ['Ltension_%s'%s for s in 'xyz']
if 0:
    fieldlist = ['Lpressure_%s'%s for s in 'x']
if 'framelist' not in dir():
    pdb.set_trace()
    framelist = range(0,110,10)
    fieldlist = ['density', 'magnetic_energy']
    fieldlist = []
    multi=True
    #fieldlist+=['velocity_magnitude']

methods={'DivB':'mip', 'Bx':'mip'}
electric = False
for sim in simlist:
    for frame in framelist:
        if extra == -1:
            name = '%s/%s/DD%04d/data%04d'%(basedir,simdict[sim],frame,frame)
            outname = 'p47_%s_n%04d'%(sim,frame)
        else:
            name = '%s/%s/ED%02d_%04d/Extra%02d_%04d'%(basedir,simdict[sim],extra,frame,extra,frame)
            outname = 'p47_%s_e%02dn%04d'%(sim,extra,frame)

        ds = yt.load(name)
        ds.periodicity=(True,True,True)
        a =( (ds.domain_right_edge-ds.domain_left_edge)/ds.domain_dimensions ).max()
        L = ds.domain_left_edge  + 3.5*a
        R = ds.domain_right_edge - 3.5*a
        C = 0.5*(L+R)
        reg = ds.region( C,L,R)
        ad=ds.all_data()
        stat(ad[('enzo','DivB')], 'divb')
        
        if multi:
            subset='_f'
            fig = plt.figure()
            grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                            nrows_ncols = (3, 2),
                            axes_pad = 1.0,
                            label_mode = "1",
                            share_all = True,
                            cbar_location="right",
                            cbar_mode="each",
                            cbar_size="3%",
                            cbar_pad="0%")
            fields=['density']
            all_new = ['Bstretching', 'Bcompression', 'Badvection', 'Ltension', 'Lpressure']

            zmindict={'Ltension_f':-100,'Lpressure_f':-100}
            zmaxdict={'Ltension_f':+100,'Lpressure_f':+100}
            fields += ["%s%s"%(f,subset) for f in all_new]
            p = yt.SlicePlot(ds, 'z', fields)
            #p.set_log('velocity_x', False)
            #p.set_log('velocity_y', False)
            for i, field in enumerate(fields):
                if field != 'density':
                    p.set_log(field,False)
                    p.annotate_magnetic_field()

                plot = p.plots[field]
                plot.figure = fig
                plot.axes = grid[i].axes
                plot.cax = grid.cbar_axes[i]
                #plot.zmin = zmindict.get(field,None)
                #plot.zmax = zmaxdict.get(field,None)
            p._setup_plots()
            this_outname = outname + "multi" + subset
            plt.savefig(this_outname)
            print this_outname
            plt.close(fig)

        for field in fieldlist:
            stat(ad[field],field)

            
            for ax in []: #'z':
                axnum = 'xyz'.index(ax)
                field_x_add = '%s'%ax
                field_y_add = '%s'%{'x':'y','y':'z','z':'x'}[ax]
                field_z_add = '%s'%{'x':'z','y':'x','z':'y'}[ax]
                #stat(reg[field], field)
                #proj= yt.ProjectionPlot(ds,ax,field, method=methods.get(field,'integrate')) #,data_source=reg)
                proj= yt.SlicePlot(ds,ax,field)
                proj.annotate_magnetic_field()
                #proj.annotate_velocity()
                #proj.annotate_streamlines('Badvection_%s'%field_y_add,'Badvection_%s'%field_z_add,factor = 16) #needs to by X and Y
                #proj.annotate_streamlines('%s-velocity'%field_y_add,'%s-velocity'%field_z_add,factor = 16) #needs to by X and Y
                #proj.annotate_streamlines('Ltension_%s'%field_y_add,'Ltension_%s'%field_z_add,factor = 16) #needs to by X and Y
                #proj.annotate_streamlines('Lpressure_%s'%field_y_add,'Lpressure_%s'%field_z_add,factor = 16) #needs to by X and Y
                #proj.set_zlim(field,-50,50)
                print proj.save(outname)
            if electric:
                dir_name = "%s/%s"%(basedir,simdict[sim])
                grid = 1
                grid_name = dsdiff_helpers.find_grid_filename(dir_name,frame,grid)
                try:
                    fptr = h5py.File(grid_name)
                    for efield in ['Ex','Ey','Ez']:
                        this_e = fptr['Grid%08d'%grid][efield][:]
                        plt.clf()
                        plot=plt.imshow( this_e.max(axis=axnum) , interpolation='nearest',origin='lower')
                        plt.colorbar(plot)
                        e_outname = "%s_%s_%s_max.png"%(outname,efield, ax)
                        plt.savefig(e_outname)
                        print e_outname
                except:
                    raise
                finally:
                    fptr.close()
