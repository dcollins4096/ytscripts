
#Plot[-D[1/2 D[sin(2 pi x)^2+cos(4 pi x)^2,y],y] -D[1/2 D[sin(2 pi x)^2+cos(4 pi x)^2,x],x] +2 D[sin( 2 pi y),y] D[ sin(4 pi x),x],{x,0,1},{y,0,1}]
if 'mhd' not in dir():
    ef('p13_fields.py')
ds = yt.load('/scratch1/dcollins/TestRunner/OT/DD0000/data0000')

proj = yt.ProjectionPlot(ds,2,'Lambda_magn')
proj.save()
