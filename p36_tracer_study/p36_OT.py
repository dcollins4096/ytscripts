
temp='/scratch1/dkl14b/orzog_tang_Mass_Interpo_bigger_Grid_Dimension/DD%04d/data%04d'
temp='/scratch1/dcollins/Paper36_tracertests/OT1/DD%04d/data%04d'
for frame in range(30):
    ds = yt.load(temp%(frame,frame))
    proj = yt.ProjectionPlot(ds,2,'density')
    proj.set_cmap('density','Greys')
    proj.annotate_streamlines('x-velocity','y-velocity',plot_args={'color':'b'})
    proj.annotate_streamlines('Bx','By',plot_args={'color':'y'})
    #proj.annotate_particles(1, col='r')
    
    print proj.save('p36c_%04d'%frame)


