


plt.clf()
map = algae_map(200)
c_to_use = [0.5]*4

alpha = at.alpha( 'x','cl2d_two','ppv_1',background_column = 2)
plot_max = 500
for n_clump, alpha_2d in enumerate(alpha):
    if at.clump_sets.has_key(n_clump):
        alpha_3d_to_use = at.clump_sets[n_clump]['Alpha2']
    else:
        continue
        alpha_3d_to_use = [0.04]
    print alpha_3d_to_use
    plot_max = max( plot_max, alpha_3d_to_use.max() )
    plt.plot([alpha_2d]*len(alpha_3d_to_use),alpha_3d_to_use,marker='o',ms=5,c=c_to_use)



if 1:
    plot_min = 0.02

    plt.plot([1,1],[plot_min,plot_max],c='k')
    plt.plot([plot_min,plot_max],[1,1],c='k')
    plt.plot([plot_min,plot_max],[plot_min,plot_max],c='k')
    plt.ylim([plot_min,plot_max])
    plt.xlim([plot_min,plot_max])
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$\alpha_{\rm{3D}}$')
    plt.xlabel(r'$\alpha_{\rm{2D}}$')
outname = 'p14b_2d_3d_%s_%04d_%s.pdf'%(at.sim,at.frame,'x')
plt.savefig(outname)
print outname
