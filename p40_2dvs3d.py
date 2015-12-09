from yt.analysis_modules.level_sets.api import * #for clumps
dir1 = '/Users/dcollins/scratch/Paper05/OK4'
dir1 = '/scratch1/dcollins/Paper05/OK4/'
ef('zeeman_measurements.py')

frame=700
size = 128*2**4
ds1 = yt.load(dir1+"/RS%04d/restart%04d"%(frame,frame))
LOS = 'x'
Blos = 'Bx'

#
# Get clumps
#
ad = ds1.all_data()
master_clump = Clump(ad,"density")
master_clump.add_validator("min_cells", 20)
c_min = ad["gas", "density"].min()
c_max = ad["gas", "density"].max()
n_steps = 5
step_size = np.exp( (np.log(c_max) - np.log(c_min))/n_steps)
find_clumps(master_clump, c_min, c_max, step_size)
bottom=get_lowest_clumps(master_clump)

#
# Plot things
#
proj = yt.ProjectionPlot(ds1,LOS,'density')
proj.set_cmap('density','gray')
proj.annotate_clumps(bottom)
print proj.save('p40_3d_test.png')


#
# Average field and density
#
#       Definitions from Collins+2011
#       Blos below is identical.  Column Density is a little different, but they should be reasonably close.  
#       for d,x in enumerate(axes):
#           tx = trans(daxes,d) #this gives [1,2] for x='x', [0,2] for x='y', [0,1] for x='z'
#           A = clump[tx[0]]*clump[tx[1]]
#           Z = clump[daxes[d]]*clump.data.convert('cm')
#           #self.AvgColumnDensity[x] = (clump['CellVolume']*clump['NumberDensity']).sum()/clump['CellVolume'].sum()
#           self.AvgColumnDensity[d] = (A*Z*clump['NumberDensity']).sum()/A.sum()
#           self.AvgBlos[d] = (clump['CellMass']*clump['MagneticField_C_%1d'%(d+1)]).sum()/clump['CellMass'].sum()


Bfield = []
ColumnDensity = []
M_list = []
R_list = []
Dx=[]
Dy=[]
Dz=[]
for n, c in enumerate(bottom):
    Dx.append( c['x'].max()-c['x'].min())
    Dy.append( c['y'].max()-c['y'].min())
    Dz.append( c['z'].max()-c['z'].min())
    R = np.sqrt(Dx[-1].v*Dy[-1].v)
    if R == 0:
        print "R=0", n
        continue
    R_list.append(R)
    A = np.pi*R**2
    M = c['cell_mass'].v.sum()
    M_list.append(M)
    Bfield.append( (c['Bz'].v*c['cell_mass'].v).sum()/M)
    ColumnDensity.append( M/A)

#
# More plotting
#
plt.clf()
plt.scatter(ColumnDensity, np.abs(Bfield))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Column Density, code')
plt.ylabel('Bfield, code')
plt.scatter(TC[1],TC[0],marker='^', c='k',label='tC2008')
plt.scatter(FT[1],FT[0],marker='D',c='c', label = 'FT2008')
plt.xlim(1e-3,1e4)
plt.ylim(0.1,1e4)
outname ='ok4_Bsigma.pdf'
plt.savefig(outname)
print outname

#
# 2d stuff
#

#ds1 = yt.load(dir1+"/RS%04d/restart%04d"%(frame,frame))
proj_flat = ds1.proj('density',LOS)
proj_weight=ds1.proj(Blos, LOS, weight_field='cell_mass')

frb_flat = proj_flat.to_frb(1,[size]*2)
frb_weight = proj_weight.to_frb(1,[size]*2)

density_flat = frb_flat['density']
cell_mass_flat = frb_flat['cell_mass']
Blos_weight = frb_weight[Blos]


data_set = {'density':the_copy(density_flat.v.reshape(size,size,1)),
            Blos+"abs":np.abs(the_copy(Blos_weight.v.reshape(size,size,1))),
            'cell_mass':the_copy(cell_mass_flat.v.reshape(size,size,1))}
bbox = np.array([[0.,1.]]*3)
ds2 = yt.load_uniform_grid(data_set, [size,size,1], length_unit="Mpc", bbox=bbox)
phase1=yt.PhasePlot(ds2.all_data(),'density',Blos+"abs",'cell_mass',weight_field=None)

#
# Plot 2d stuff
#
outname = 'p40_badone_abs.pdf'
phase1.save(outname)
print outname

#
# Plot cores and projections together.
#

ploot = phase1.plots[('gas','cell_mass')].axes
ploot.scatter(ColumnDensity,np.abs(Bfield),c='r')
ploot.scatter(TC[1],TC[0],marker='^', c='k',label='tC2008')
ploot.scatter(FT[1],FT[0],marker='D',c='c', label = 'FT2008')
outname = 'p40_2d3d.pdf'
phase1.save(outname)
print outname
