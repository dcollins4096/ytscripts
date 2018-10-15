
if 'g' not in dir():
    g= ds.index.grids[0]

d  = g[('enzo','Density')]
te = g[('enzo','TotalEnergy')]*d
vx = g[('enzo','x-velocity')]
vy = g[('enzo','y-velocity')]
vz = g[('enzo','z-velocity')]
ke = 0.5*d*(vx*vx+vy*vy+vz*vz)
bx = g[('enzo','Bx')]
by = g[('enzo','By')]
bz = g[('enzo','Bz')]
be = 0.5*(bx*bx+by*by+bz*bz)
