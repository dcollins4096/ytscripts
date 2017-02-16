if 0:
    g,g2 = fltc('output.append(car.ds.index.grids[-1])')
    II = (g['x']-g.LeftEdge[0]-0.5*g.dds[0])/g.dds[0]
    JJ = (g['y']-g.LeftEdge[1]-0.5*g.dds[1])/g.dds[1]
    KK = (g['z']-g.LeftEdge[2]-0.5*g.dds[2])/g.dds[2]
    index = (II + g.shape[0]*(JJ+g.shape[1]*KK)).v
    II = (g2['x']-g2.LeftEdge[0]-0.5*g2.dds[0])/g2.dds[0]
    JJ = (g2['y']-g2.LeftEdge[1]-0.5*g2.dds[1])/g2.dds[1]
    KK = (g2['z']-g2.LeftEdge[2]-0.5*g2.dds[2])/g2.dds[2]
    index2= (II + g2.shape[0]*(JJ+g2.shape[1]*KK)).v

metalfields = [
   "Cooling_Time",
   "DebugField",
   "Density",
   "Electron_Density",
   "GasEnergy",
   "HII_Density",
   "HI_Density",
   "HeIII_Density",
   "HeII_Density",
   "HeI_Density",
   "Metal_Density",
   "Temperature",
   "TotalEnergy",
   "x-velocity",
   "y-velocity",
   "z-velocity",'kinetic_energy']
fields=['mjeans','Temperature','jeans_density','density']+metalfields
mask = g['bmass']>g['mjeans']
for field in fields:
    plt.clf()
    print field
    f1 = g[field][mask].flatten()
    f2 = g2[field][mask].flatten()
    value = (f2-f1)/(0.5*(f1+f2))
    colors = nar([[0,1,0,0.5]]*len(value))
    colors[value < 0] = [1,0,0,0.5]
    value = np.abs(value)
    plt.scatter(f1,f2,c=colors,linewidth=0)
    plt.xlabel(field + ' ppm')
    #plt.ylabel('difference')
    plt.ylabel(field + 'hydro6')
    plt.xscale('log')
    plt.yscale('log')
    mm=min([f1.min(),f2.min()])
    xx=max([f1.max(),f2.max()])
    plt.xlim(mm,xx)
    #plt.ylim(value.min(), value.max())
    plt.ylim(mm,xx)
    plt.plot([mm,xx],[mm,xx])
    plt.savefig('p33_dumbplotb_%s.pdf'%field)
