import matplotlib.colors as colors
plt.close('all')

for frame in [0,20,50,80,100]:
    for car in flt.taxi_list:
        ds=car.load(frame)
        level=0
        cg = ds.covering_grid(level,ds.domain_left_edge, ds.domain_dimensions)
        plt.clf()
        Bx=cg['Bx'].v
        By=cg['By'].v
        Bz=cg['Bz'].v
        Bx0 = Bx.mean()
        By0 = By.mean()
        Bz0 = Bz.mean()
        if 1:
            BxF = Bx-Bx0
            print(BxF.mean(), "xxxxx")
            plt.hist(2*Bx.flatten(), color='r',label='bx',histtype='step')
            plt.hist(BxF.flatten(), color='k',label='bx -bx0',histtype='step')
            plt.hist(0.5*By.flatten(), color='g',label='by',histtype='step')
            plt.hist(Bz.flatten(), color='b',label='bz',histtype='step')
            #plt.hist(np.log10(np.abs(Bx.flatten())), color='r',label='bx',histtype='step')
            #plt.hist(np.log10(np.abs(BxF.flatten())), color='k',label='bx -bx0',histtype='step')
            #plt.hist(np.log10(np.abs(By.flatten())), color='g',label='by',histtype='step')
            #plt.hist(np.log10(np.abs(Bz.flatten())), color='b',label='bz',histtype='step')
            plt.legend(loc=0)

        plt.yscale('log')
        plt.savefig('Bhist.png')
        plt.clf()
        Bxflat=(Bx.flatten())[::10]
        Byflat=(By.flatten())[::10]
        ok = np.logical_and(Bxflat != 0, Byflat != 0)


        fig,ax=plt.subplots(1,1, figsize=(5.5,5))
        ax.set_aspect('equal')
        axis=1
        TheX = cg['z'].mean(axis=axis).v/1e5
        TheY = cg['x'].mean(axis=axis).v/1e5

        TheBx = Bz-Bz0
        TheBy = Bx-Bx0
        TheBxProj = (TheBx).sum(axis=axis) 
        TheByProj = (TheBy).sum(axis=axis)
        R = np.sqrt(TheX**2+TheY**2)



        Ni =  cg['Density_56Ni'].v.max(axis=axis)
        if (Ni == 0).any():
            MinMin = Ni[ Ni>0].min()
        else:
            MinMin = Ni.min()
        MinMin = 1e6
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
        plt.savefig('Bquiv_%s_n%04d.png'%(car.name,frame))

    if 0:

        fig,ax=plt.subplots(2,2)
        axlist=ax.flatten()
        for  aaa in axlist:
            aaa.set_aspect('equal')
        axis=1
        TheX = cg['z'].mean(axis=axis).v/1e5
        TheY = cg['x'].mean(axis=axis).v/1e5
        axis=1
        TheX = cg['z'].mean(axis=axis).v/1e5
        TheY = cg['x'].mean(axis=axis).v/1e5

        TheBx = Bz#-2*Bz0
        TheBy = Bx#-2*Bx0
        TheBxProj = (TheBx).sum(axis=axis) 
        TheByProj = (TheBy).sum(axis=axis)
        norm = colors.SymLogNorm(vmin=-MaxMax,vmax=MaxMax, linthresh = 1e4)
        #cmap = copy.copy(mpl.cm.get_cmap("viridis"))
        px=ax[0][0].pcolormesh(TheX, TheY, TheBxProj)#,norm=norm,cmap='jet')
        py=ax[1][0].pcolormesh(TheX, TheY, (np.abs(TheByProj)))#,norm=norm,cmap='jet')
        ptheta=ax[0][1].pcolormesh(TheX, TheY, np.arctan2(TheByProj,TheBxProj),cmap='twilight')
        fig.colorbar(px,ax=ax[0][0])
        fig.colorbar(py,ax=ax[1][0])
        fig.colorbar(ptheta,ax=ax[0][1])
        R = np.sqrt(TheX**2+TheY**2)
        ax[1][1].streamplot(np.unique(TheX), np.unique(TheY), TheBxProj, TheByProj ,linewidth=0.1,color='b', arrowsize=0.1)
        fig.savefig('Btest.png')
