
if 1:
    bottom=get_lowest_clumps(master_clump)
    proj = yt.ProjectionPlot(ds1,0,'density')
    proj.set_cmap('density','gray')
    proj.annotate_clumps(bottom)
    print proj.save('p40_3d_test.png')

if 1:
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
ef('zeeman_measurements.py')
if 1:
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


