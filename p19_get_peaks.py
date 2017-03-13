from yt.analysis_modules.level_sets.api import * #for clumps
from yt.utilities.data_point_utilities import FindBindingEnergy

frame = 95; basedir = '/scratch1/dcollins/Paper08/B02/256/'
frame = 60; basedir = '/scratch1/dcollins/Paper08/B02/512/'
frame = 60; basedir = '/scratch1/dcollins/Paper08/B20/512/'
frame = 68; basedir = '/scratch1/dcollins/Paper08/B20/512/'; setname = 'b20'
fname = '%s/RS%04d/restart%04d'%(basedir, frame,frame)

if 1:
  ds = yt.load(fname)
  region = ds.region([0.5]*3, [0.25]*3, [0.75]*3)
  master_clump = Clump(region,"density")
  t0 = time.time()
  find_clumps(master_clump, 5e5,1e8,10)
  t1 = time.time()
  print "took", t1-t0 # about 10 minutes.

if 1:
    leaf_clumps = get_lowest_clumps(master_clump)
    proj = yt.ProjectionPlot(ds,0,'density')
    proj.set_cmap('density','gray')
    proj.annotate_clumps(leaf_clumps)
    proj.save('p24_leaves_t1.png')

if 1:
    peak_list = []
    for cl in leaf_clumps:
        peak_ind = np.where( cl['density'] == cl['density'].max())[0][0]
        xpeak = cl['x'][peak_ind]
        ypeak = cl['y'][peak_ind]
        zpeak = cl['z'][peak_ind]
        peak_list.append( nar([xpeak,ypeak,zpeak]))

    proj = yt.ProjectionPlot(ds,0,'density')
    for n,peak in enumerate(peak_list):
        proj.annotate_point(peak,n,text_args={'color':'r'})
    proj.save('b025v_n0060_peaks_2.png')
    fPickle.dump(peak_list,"p12_%s_512_n%04d_peaks2.pickle"%(setname,frame))



