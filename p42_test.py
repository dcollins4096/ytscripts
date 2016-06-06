if  'ef' not in dir():
    execfile('go')
basedir = '/Users/dcollins/scratch/Paper42_NewForcing/r01_test'
if 'frame' not in dir():
    frame = 0
for frame in [0,5]:
    ds = yt.load('%s/DD%04d/data%04d'%(basedir,frame,frame))
    g0 = ds.index.grids[0]
    d= g0['density'].v
    db=g0['DivB'].v
    stat(db,frame)
    #v0 = v[0,0,0]
    #stat(v-v[0], frame)
    #print "v0", v0
#frame+=1
#print yt.ProjectionPlot(ds,0,'x-velocity').save()
