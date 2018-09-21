
from yt.analysis_modules.level_sets.api import * #for clumps
frame = 125
scratchdir = '/scratch1/dcollins/Paper19/u05-r4-l4-128'
scratchdir = '/work/00369/tg456484/maverick/Paper19/u05-r4-l4-128'
fname = '%s/DD%04d/data%04d'%(scratchdir,frame,frame)

if 'ds' not in dir():
    ds = yt.load(fname)

if 'loc' not in dir():
    val, loc = ds.find_max('density')

if 0:
    locn = loc.to_ndarray()
    Left = nar([locn[0]-0.3, 0, locn[2]+0.2])
    Right = nar([locn[0]+0.1, 1, locn[2]+0.5])
    center = 0.5*(Left+Right)
    rect = ds.region(center,Left,Right)
    for ax in 'xyz':
        proj = ds.proj('density', ax, data_source = rect, center = loc)
        pw = proj.to_pw(center = loc)
        pw.annotate_grids()
        print pw.save('u05_0125_region_withgrids')

if 0:
    for ax in 'xyz':
        print "proj", ax
        proj = yt.ProjectionPlot(ds,ax,'density', center=loc)
        #proj.annotate_grids()
        print proj.save('u05_0125_nogrids')

if 0:
    prof = yt.ProfilePlot(ds.all_data(),'density','cell_volume',weight_field=None)
    prof.save('u05_0125')

if 0:
    master_clump = Clump(ds.all_data(),"density")
    master_clump.add_validator("min_cells", 20)
    c_min = 10
    c_max = 1e7
    step = 10
    find_clumps(master_clump, c_min, c_max, step)

if 0:
    if 'leaf_clumps' not in dir():
        leaf_clumps = get_lowest_clumps(master_clump)
    for ax in 'xyz':
        print "proj", ax
        proj = yt.ProjectionPlot(ds,ax,'density')
        proj.annotate_clumps(leaf_clumps)
        print proj.save('u05_0125_leaf_clumps')

if 0:
    for c in leaf_clumps:
        ind= c['density'].max()
        
if 0:
    width = ds.arr(0.05,'code_length')
    if 'loc' not in dir():
        val,loc = ds.find_max('density')
    sphere = ds.sphere(loc,width)
    other_master_clump = Clump(sphere,"density")
    other_master_clump.add_validator("min_cells", 20)
    c_min = 10
    c_max = 1e7
    step = 10
    find_clumps(other_master_clump, c_min, c_max, step)



if 0:
    if 'other_master_clump' not in dir():
        other_master_clump = master_clump
    leaf_clumps = get_lowest_clumps(other_master_clump)
    print "arf"

    def _peak(clump):
        print "c",
        peak_ind = np.where( clump['density'] == clump['density'].max())[0][0]
        xpeak = clump['x'][peak_ind]
        ypeak = clump['y'][peak_ind]
        zpeak = clump['z'][peak_ind]
        xmin = clump['x'].min(); xmax=clump['x'].max()
        ymin = clump['y'].min(); ymax=clump['y'].max()
        zmin = clump['z'].min(); zmax=clump['z'].max()
        return "Peak (%.6e, %.6e, %.6e)  x (%.6e, %.6e) y (%.6e, %.6e) z (%.6e, %.6e)"%(xpeak, ypeak, zpeak, xmin, xmax, ymin, ymax, zmin, zmax)
    add_clump_info("peak", _peak)
    other_master_clump.add_info_item("peak")
    write_clumps(other_master_clump,0, "%s_smallclumps.txt" % ds)

if 0:
    peak_list = []
    for cl in leaf_clumps:
        peak_ind = np.where( cl['density'] == cl['density'].max())[0][0]
        xpeak = cl['x'][peak_ind]
        ypeak = cl['y'][peak_ind]
        zpeak = cl['z'][peak_ind]
        peak_list.append( nar([xpeak,ypeak,zpeak]))
        #xmin = cl['x'].min(); xmax=cl['x'].max()
        #ymin = cl['y'].min(); ymax=cl['y'].max()
        #zmin = cl['z'].min(); zmax=cl['z'].max()

if 0:
  master_clump = Clump(ds.all_data(),"density")
  master_clump.add_validator("min_cells", 20)
  c_min = 10
  c_max = 1e7
  step = 10
  find_clumps(master_clump, c_min, c_max, step)

if 1:
    if 'peak_list' not in dir():
        peak_list = fPickle.load('u05_0125_peaklist.pickle')
    keepers = [0,1,8,10,11,12,67,64,61, 201, 125, 306]
    for ax in 'xyz':
        print "proj", ax
        proj = yt.ProjectionPlot(ds,ax,'density')
        proj.set_cmap('density','gray')
        #for n, peak in enumerate(peak_list):
        for n in keepers:
            peak = peak_list[n]
            proj.annotate_point(peak,n, text_args={'color':'r'})
        proj.save('u05_0125_all_peaks_keepers')
        #proj.annotate_grids()
        

