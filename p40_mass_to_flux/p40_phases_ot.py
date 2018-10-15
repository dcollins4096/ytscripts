
proj = ds1.proj('density',2)
frb = proj.to_frb(1,16)
proj2 = ds1.proj('Bx',2,weight_field='cell_mass')
frb2 = proj2.to_frb(1,16)
density = frb['density'].v.reshape(16,16,1)
data = {'density':density, 'Bx':frb2['Bx'].reshape(16,16,1), 'cell_mass':frb['cell_mass'].reshape(16,16,1)}
bbox = np.array([[0.,1.]]*3)
ds = yt.load_uniform_grid(data, [16,16,1], length_unit="Mpc", bbox=bbox, nprocs=64)
ad=ds.all_data()
plot = yt.PhasePlot(ad,'density','Bx','cell_mass',weight_field=None)
plot.save('test.png')
