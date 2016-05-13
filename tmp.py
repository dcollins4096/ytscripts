if 0:
    def printer(at1,field):
        for n in sorted(at.clump_sets[3][field]):
            print "%0.2e"%n

    thing2 = at.clump_sets[3]
    M  = thing2['Mass']
    v2 = thing2['VelocityDispersion']
    r2 = thing2['R2']
    G = 1620./(4*np.pi)
    a1 = 5.*v2*v2*r2/(3*M*G)
    a2 = thing2['Alpha2']
    for n,m in enumerate(thing2['Mass']):
        print " m %0.2e v %0.2e r %0.2e a1 %0.2e a2 %0.2e"%(m,v2[n],r2[n],a1[n],a2[n])
if 1:
    for axis in [0,1,2]:
        proj = yt.ProjectionPlot(at.ds,axis,'Density')
        proj.annotate_clumps(at.leaf_clumps[3])
        print proj.save('p14b_%s_yt3_n%04d_cl2d_three_cl%02d'%(at.sim,at.frame,3))
