if 'ef' not in dir():
    execfile('go')
import dt_hierarchy
reload(dt_hierarchy)

flt = taxi.fleet(['aj21','aj22'])
fname1, fname2 = flt('output += [car.directory+"/dump"]')
level=4
if 1:
#fname1= '/scratch1/dcollins/Paper33_galaxies/aj19_r2_ppm/dump'; f1=aj19
#fname2= '/scratch1/dcollins/Paper33_galaxies/aj20_r2_h6/dump'; f2=aj19
    vf1=dt_hierarchy.vis_file(fname1)
    vf2=dt_hierarchy.vis_file(fname2)
    prefix = flt.allnames()

    plt.clf()
    plt.plot(vf1.time[vf1.l[level]], vf1.dt[vf1.l[level]],c='r',label='aj19')
    plt.plot(vf2.time[vf2.l[level]], vf2.dt[vf2.l[level]],c='g',label='aj20')
#plt.hist(vf1.dt[vf1.l[4]],histtype='step')
#plt.hist(vf2.dt[vf2.l[4]],histtype='step')
    outname = '%sl%s_dttime.png'%(prefix,level)
    plt.savefig(outname)
    print outname
    plt.clf()
    plt.hist(vf1.dt[vf1.l[4]],histtype='step')
    plt.hist(vf2.dt[vf2.l[4]],histtype='step')
    outname = '%s_l%s_hist.png'%(prefix,level)
    plt.savefig(outname)
    print outname
    if 0:
        mslice = slice(None)
        plt.plot(TheX,TheY ,c=[0.5]*3)
        plt.scatter(TheX, TheY,marker='*',c=vf1.color_list,linewidths=0, s=100)
        plt.xlabel(TheXLabel); plt.ylabel(TheYLabel);
        plt.yscale('log')
        plt.savefig(output)
        print output
