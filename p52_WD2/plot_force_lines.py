import matplotlib.colors as colors
plt.close('all')

for frame in [0]: #[0,20,50,80,100]:
    for car in flt.taxi_list[3:4]:
        ds=car.load(frame)
        ad=ds.all_data()

        level=0
        cg = ds.covering_grid(level,ds.domain_left_edge, ds.domain_dimensions,num_ghost_zones=1)
        plt.clf()

        field = 'Lpressure'
        #field = 'Ltension'
        #field = 'magnetic_field'
        Bx=cg['%s_x'%field].v
        By=cg['%s_y'%field].v
        Bz=cg['%s_z'%field].v
        Bx0 = Bx.mean()
        #By0 = By.mean()
        Bz0 = Bz.mean()

        fig,ax=plt.subplots(1,1, figsize=(5.5,5))
        ax.set_aspect('equal')
        axis=1
        TheX = cg['z'].mean(axis=axis).v/1e5
        TheY = cg['x'].mean(axis=axis).v/1e5
        R = np.sqrt(TheX**2+TheY**2)

        TheBx = Bz#-Bz0
        TheBy = Bx#-Bx0
        TheBxProj = (TheBx).sum(axis=axis) 
        TheByProj = (TheBy).sum(axis=axis)
        #Ny = int(TheBx.shape[1]//2)
        #TheBxProj = TheBx[:,Ny,:]
        #TheByProj = TheBy[:,Ny,:]

        Ni =  cg['Density_56Ni'].v.max(axis=axis)
        if (Ni == 0).any():
            MinMin = Ni[ Ni>0].min()
        else:
            MinMin = Ni.min()
        MaxMax = Ni.max()
        norm = colors.LogNorm(vmin=MinMin,vmax=MaxMax)
        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')
        plot=ax.pcolormesh(TheX,TheY, Ni,norm=norm,cmap=cmap)

        Rho =  (cg['Density'].sum(axis=axis).v/cg['Density'].shape[axis])
        Rho = cg['Density'].max(axis=axis)
        levels = [9.3e7, 5e8, 1e9]#, 2e8,9e8,1e9]

        ax.set_title(r"$t=%0.1f\ \rm{s}$"%(frame/100))
        ax.streamplot(np.unique(TheX), np.unique(TheY), TheBxProj, TheByProj ,linewidth=0.1,color='b', arrowsize=0.1,density=2)

        circle = plt.Circle( (0,0), 1975, edgecolor='k', facecolor=[0,0,0,0])
        ax.add_patch(circle)
        ax.set_xlabel(r'$z\ [\rm{km}]$')
        ax.set_ylabel(r'$x\ [\rm{km}]$')
        fig.savefig('Fquiv_%s_%s_n%04d.png'%(field,car.name,frame))


        fig2,ax2=plt.subplots(1,1, figsize=(5.5,5))
        #field='magnetic_pressure'

        #field = 'Ltension_mag'
        #tots = cg[field]
        tots = np.abs(Bz)#Bx**2+By**2#+Bz**2
        array=TheByProj**2+TheBxProj**2 #tots.sum(axis=axis)
        #array=tots.sum(axis=axis)
        if 0:
            MinMin = tots[ tots>0].min()
            MaxMax = tots.max()/1e3
        if 0:
            MinMin = 1e14
            MaxMax = 1e37
            ax2.contour(TheX,TheY,array,levels=[1e13,1e14,1e15, 1e16, 1e17, 1e18])
        if 0:
            MinMin = 1e14
            MaxMax = 1e37
            ax2.contour(TheX,TheY,array,levels=[1e30, 1e31,1e32,1e35])
        if 1:
            MinMin = array[ array>0].min()
            MaxMax = array.max()
            ax2.contour(TheX,TheY,array,levels=[1e30, 1e31,1e32,1e35])
            #ax2.streamplot(TheX[0,:], TheY[:,0], TheBxProj, TheByProj ,linewidth=0.1,color='b', arrowsize=0.1,density=2)
            ax2.streamplot(TheX, TheY, TheBxProj, TheByProj ,linewidth=0.1,color='b', arrowsize=0.1,density=2)

        if 0:
            MinMin = array[ array>0].min()
            MaxMax = array.max()

        norm = colors.LogNorm(vmin=MinMin,vmax=MaxMax)
        #cmap = copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')
        plot=ax2.pcolormesh(TheX,TheY, array,norm=norm,cmap=cmap)
        fig2.colorbar(plot,ax=ax2)
        fig2.savefig("Fcmb_%s_%s_n%04d.png"%(field,car.name,frame))



