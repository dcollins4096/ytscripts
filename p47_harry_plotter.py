if 'ef' not in dir():
    execfile('go')

import dsdiff_helpers
import dsdiff
reload(dsdiff)
basedir = '/Users/dcollins/scratch/Paper47_forces/'
simdict={}
simdict['a01'] = 'a01'
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
if 1:
    fieldlist = ['Lpressure_%s'%s for s in 'x']

methods={'DivB':'mip', 'Bx':'mip'}
electric = False
for sim in ['a02']:
    for frame in framelist:
        if extra == -1:
            name = '%s/%s/DD%04d/data%04d'%(basedir,sim,frame,frame)
            outname = '%s_n%04d'%(sim,frame)
        else:
            name = '%s/%s/ED%02d_%04d/Extra%02d_%04d'%(basedir,sim,extra,frame,extra,frame)
            outname = '%s_e%02dn%04d'%(sim,extra,frame)

        ds = yt.load(name)
        ds.periodicity=(True,True,True)
        a =( (ds.domain_right_edge-ds.domain_left_edge)/ds.domain_dimensions ).max()
        L = ds.domain_left_edge  + 3.5*a
        R = ds.domain_right_edge - 3.5*a
        C = 0.5*(L+R)
        reg = ds.region( C,L,R)
        ad=ds.all_data()
        stat(ad[('enzo','DivB')], 'divb')

        for field in fieldlist:
            stat(ad[field],field)
            
            for ax in 'x':
                axnum = 'xyz'.index(ax)
                field_x_add = '%s'%ax
                field_y_add = '%s'%{'x':'y','y':'z','z':'x'}[ax]
                field_z_add = '%s'%{'x':'z','y':'x','z':'y'}[ax]
                #stat(reg[field], field)
                #proj= yt.ProjectionPlot(ds,ax,field, method=methods.get(field,'integrate')) #,data_source=reg)
                proj= yt.SlicePlot(ds,ax,field)
                #proj.annotate_magnetic_field()
                #proj.annotate_streamlines('Badvection_%s'%field_y_add,'Badvection_%s'%field_z_add,factor = 16) #needs to by X and Y
                #proj.annotate_streamlines('%s-velocity'%field_y_add,'%s-velocity'%field_z_add,factor = 16) #needs to by X and Y
                #proj.annotate_streamlines('Ltension_%s'%field_y_add,'Ltension_%s'%field_z_add,factor = 16) #needs to by X and Y
                proj.annotate_streamlines('Lpressure_%s'%field_y_add,'Lpressure_%s'%field_z_add,factor = 16) #needs to by X and Y
                #proj.set_zlim(field,-50,50)
                print proj.save(outname)
            if electric:
                dir_name = "%s/%s"%(basedir,sim)
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
