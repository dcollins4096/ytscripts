if 'ef' not in dir():
    execfile('go')
import dsdiff_helpers
import dsdiff
reload(dsdiff)

#basedir = '/Users/dcollins/scratch/EnzoProjects/E21_AdiabaticCorners'
basedir = '/Users/dcollins/scratch/EnzoProjects/E12/'
e01 = 'e01_ae_ppm'
e02 = 'e02_ae_ded'
e03 = 'e03_ae_ct'
e04 = 'e04_ae_ded_old'
e05 = 'e05_ae_amr_ppm'
e06 = 'e06_ae_amr_ded'
e07 = 'e07_ae_amr_ct'
e12 = 'e12_PR_ded_amr'

sim1 = e12
#sim1 = 'a08_post_sbc'
if 'frame' not in dir():
    frame = 16

dir_name1 = '%s/%s'%(basedir,sim1)
dir_name2 = '%s/%s'%(basedir,sim2)
fields = adiabatic_hydro
if 'field' not in dir():
    field = 'Density'

if 0:
    if 'grid' not in dir():
        grid = 1
    fname = dsdiff_helpers.find_grid_filename(dir_name1,frame,grid)
    print fname
    fptr = h5py.File(fname)
    try:
        data = fptr["Grid%08d"%grid][field][:]
        mean = data[0,0,0]
        if mean == 0:
            over = 1
        else:
            over=mean
        stat((data-mean)/over)
    except:
        raise
    finally:
        fptr.close()


if 1:
    ds_name = dsdiff_helpers.get_ds_name(dir_name1,frame)
    print ds_name
    ds = yt.load(ds_name)
    ad = ds.all_data()
    #stat(1.-ad['density'].in_units('code_density').v, sim1)
    proj = yt.ProjectionPlot(ds,0,'TotalEnergy')
    proj.annotate_grids()
#proj.set_zlim('density',zlim[0],zlim[1])
    print proj.save('%s_%04d'%(sim1,frame))
    slice=yt.SlicePlot(ds,0,('enzo','TotalEnergy'), center=[0.6,0.5,0.5])
    slice.annotate_grids()
    print slice.save('%s_%04d'%(sim1,frame))
