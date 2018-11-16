
exec(open('p52_WD2/potentials.py').read(),globals(),locals())
pot = potential(a01)
plt.clf()
plt.close('all')
if 1:
    for frame in range(1,6):
        ds = a01.load(frame)
        ad=ds.all_data()
        sph = ds.sphere([0.5]*3,0.25)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.scatter(sph['radius'],sph['density'])
        ax2.set_xlabel('r')
        ax2.set_ylabel('d')
        outname = 'p52_density_%04d.png'%frame
        fig2.savefig(outname)
        print(outname)
if 1:
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    x0 = pot.read(frame=1)
    for frame in range(1,6):
        x=pot.read(frame=frame)
        ax1.plot(x.Rb,x0.Mi-x.Mi,label='M<',marker='*')
        ax1.plot(x.Rb,x0.Ms-x.Ms,label='dM',marker='*')
        #ax1.legend(loc=2)

    outname = 'p52_test_multi_b.png'
    fig1.savefig(outname)
    print(outname)
