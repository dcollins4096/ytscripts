
f1 = '/Users/dcollins/scratch/Paper33_galaxies/e04_64_ngolbaum'
f2 = '/Users/dcollins/scratch/Paper33_galaxies/e02_64'

d1 = yt.load('%s/DD%04d/DD%04d'%(f1,38,38))
d2 = yt.load('%s/DD%04d/DD%04d'%(f2,38,38))

for ax in [0,1,2]:
    proj = yt.ProjectionPlot(d1,ax,'density')
    proj.annotate_grids()
    print proj.save('e04')
    proj=yt.ProjectionPlot(d2,ax,'density')
    proj.annotate_grids()
    print proj.save('e02')

