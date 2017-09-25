if 'ef' not in dir():
    execfile('go')
import ytree

reload(taxi)
del car
if 'car' not in dir():
    reload(taxi)
    car=taxi.taxi('c05')

if 'mytree' not in dir():
    tree_file = car.directory+"/rockstar_halos/trees/tree_0_0_0.dat"
    mytree=ytree.load(tree_file)
    plt.hist(np.log10(mytree['mvir']))
    plt.ylabel('N')
    plt.xlabel(r'$\log_{10} M_{vir}$')
    outname = '%s_z=0_halo.pdf'%car.outname
    plt.savefig(outname)


if 1:
    big_id = np.where(mytree['mvir'] == mytree['mvir'].max())[0][0]
    big_tree = mytree[big_id]


if 1:
    car.cmap['density']='gray'
    car.operation == 'Full'
    car.outname = car.name+"_halotest5"
    car.callbacks=['z_title','halos']
    #car.callback_args['spheres']={'centers':cen,'radii':rad,'ids':ids, 'circle_args':{'color':'yellow'}}
    car.name_syntax='outputlog'
    #car.frames='last'
    #print car.plot()


    last_ok=None
    big_tree = mytree
    ok = {}
    last_ok = {}
    id_offset = min(halo['tree_id'] for halo in mytree)
    for frame in [56]:
        car.frames=[frame]
        car.fill()
        for n,big_tree in enumerate(mytree):
            ok[n] = nar(np.abs(big_tree['tree','redshift'] -car.ds['CosmologyCurrentRedshift'])) < 1e-6
            if ok[n].sum() > 0 or not last_ok.has_key(n):
                these_halos = nar(big_tree['tree'])[ok[n]]
                last_ok[n]= ok[n]
            if n == 0:
                good_halos = []
            good_halos = good_halos + list(these_halos)
            #print "many", n, len(good_halos)

        car.callback_args['halos']={'halos':good_halos, 'circle_args':{'color':'yellow'}, 'id_offset':id_offset}
#       car.cmap['density']='gray'
        car.plot()

