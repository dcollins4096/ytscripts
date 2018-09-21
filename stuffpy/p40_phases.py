dir1 = '/Users/dcollins/scratch/Paper05/OK4'
dir1 = '/scratch1/dcollins/Paper05/OK4/'
frame=700
size = 128*2**4
ds1 = yt.load(dir1+"/RS%04d/restart%04d"%(frame,frame))
proj = ds1.proj('density',2)
frb = proj.to_frb(1,size)
proj2 = ds1.proj('Bx',2,weight_field='cell_mass')
frb2 = proj2.to_frb(1,size)
data = {'density':copy.copy(frb['density'].v).reshape(size,size,1),
        'Bx':copy.copy(frb2['Bx'].v).reshape(size,size,1),
        'cell_mass':copy.copy(frb['cell_mass'].v).reshape(size,size,1)}
bbox = np.array([[0.,1.]]*3)

ds = yt.load_uniform_grid(data, [size,size,1], length_unit="Mpc", bbox=bbox)
ad=ds.all_data()
plot = yt.PhasePlot(ad,'density','Bx','cell_mass',weight_field=None)
plot.save('p40_uniform_phase_different_weights.png')

#ds2 = frb.export_dataset(fields=['density','Bx','cell_mass'])
#ad2 = ds2.all_data()
#plot2 = yt.PhasePlot(ad2,'density','Bx','cell_mass',weight_field=None)
#plot2.save('p40_exported_phase.png')
#plot3 = yt.PhasePlot(ad2,'density','Bx','cell_mass',weight_field=None)
#plot2.save('p40_exported_phase.png')


