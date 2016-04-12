
if 0:
    if 'frame' not in dir():
        frame = 28
    basedir = '/scratch1/dcollins/Paper20/a05_OT'
    basedir = '/scratch1/dcollins/Paper20/a06_weight_test'; frame = 0
    prefix = 'a06p_%04d'%frame
    setname = '%s/DD%04d/data%04d'%(basedir,frame,frame)
    line_of_sight = 2
    ds = yt.load(setname)
    proj= yt.ProjectionPlot(ds,line_of_sight,'density')
    print proj.save(prefix)
    proj_mag = ds.proj('Bx',line_of_sight, weight_field = 'cell_mass') 
    print proj_mag.to_pw().save(prefix+"mass_weight")
    proj_mag_noweight = ds.proj('Bx',line_of_sight)
    print proj_mag_noweight.to_pw().save(prefix+"no_weight")

if 1:
    basedir = '/scratch1/dcollins/Paper08/B02/512'
    if 'frame' not in dir():
        frame = 30
    setname = basedir + "/RS%04d/restart%04d"%(frame,frame)
    prefix = 'B02_512_%04d'%frame
    ds = yt.load(setname)
    proj= yt.ProjectionPlot(ds,line_of_sight,'density')
    print proj.save(prefix)
    proj_mag = ds.proj('Bx',line_of_sight, weight_field = 'cell_mass')
    print proj_mag.to_pw().save(prefix+"mass_weight")
    proj_mag_noweight = ds.proj('Bx',line_of_sight)
    print proj_mag_noweight.to_pw().save(prefix+"no_weight")
