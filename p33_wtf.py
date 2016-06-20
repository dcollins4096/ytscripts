
d1 = '/Users/dcollins/scratch/Paper33_galaxies/e01_stock'

for frame in range(241):
    fname = '%s/DD%04d/DD%04d'%(d1,frame,frame)

    ds = yt.load(fname)
    for ax in [0,1,2]:
        print yt.ProjectionPlot(ds,ax,'density').save()
