
import matplotlib.colors as colors
basedir = '/home/bhristov/mytest/8enzo-dev-bh'
basedir = '/home/dcollins4096/ASTRO_12_3/runs/'
basedir = '/scratch1/dcollins/Paper10/'

frames = {182:600,181:600, 150:600,190:600}
Qi = 'QInstantaneous'

#def ni_frac(field,data):
#    return data['Density_56Ni']/data['density'].in_units('code_density').v
#yt.add_field('56NiMassFrac', function=ni_frac, units = "", take_log=False)
#field  = '56NiMassFrac'
field = 'density'
if 1:
    ds={}
    cg={}
    plt.close('all')
    fig, ax = plt.subplots()
    for run in  [190]: #[181,182,150]:
        frame=frames[run]
        ds[run] = yt.load('%s/%s/DD%04d/data%04d'%(basedir,run,frame,frame))
        size = ds[run].domain_right_edge-ds[run].domain_left_edge
        dims = ds[run].domain_dimensions
        cg[run] = ds[run].covering_grid(0,[0.0,0.0,0.0],dims)
        axis=1
        array = np.log(np.sum(cg[run][field].v,axis=axis))
        #norm = norm=colors.LogNorm(vmin=array.min(), vmax=array.max()),
        #norm = norm=colors.Normalize(vmin=array.min(), vmax=array.max())
        cax = ax.imshow(np.log(array), interpolation='nearest') #,cmap='Greys') #, norm=norm)
        #fig.colorbar(cax)

        proj = yt.ProjectionPlot(ds[run],axis,field)
        #proj.set_cmap('QInstantaneous','Greys')
        print proj.save('p10p_r%d_n%d.pdf'%(run,frame))
        outname = 'p10p_r%d_n%04d_%s_%s.pdf'%(run,frame,field,'xyz'[axis])
        fig.savefig(outname)
        print outname


if 0:
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    run = 181
    ax1.imshow(np.log(np.sum(cg[run][Qi],axis=1)), cmap='Greys')
    ax2 = fig.add_subplot(2,1,2)
    run = 182
    ax2.imshow(np.log(np.sum(cg[run][Qi],axis=1)), cmap='Greys')
    fig.savefig('test2.png')
