execfile('p42_pg.py')
import fourier_tools.fourier_filter as Filter


#aq32 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq32_m2.9_drive1_32',name='aq32',frames=range(0,42,2),fields=['density'])
#aq31 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq31_m2.9_drive0.5_32',name='aq31',frames=range(0,42,2),fields=['density'])
#aq16 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq16_m2.9_drive0_noamr_128',name='aq15',frames=range(17),fields=['density'])
#aq42 = taxi.taxi(directory='/scratch/00369/tg456484/Paper42_NewAK/aq42_m9_drive0.5_512_p59_ppm',name='aq42',frames=range(90),fields=['density'])
#aw11 = taxi.taxi('aw11')
#car = aw11
#frame = 120
#car.fill(frame)
#ds=car.ds

#fields = ['%s-velocity'%s for s in 'xyz']
#fields = ['DrivingField%s'%s for s in '123']
def all_the_spectra(car,fields, frames=None):
    if frames is None:
        frames = car.frames
    for frame in frames:
        car.fill(frame)
        cube = car.ds.covering_grid(0, left_edge=car.ds.domain_left_edge,
                                dims=car.ds.domain_dimensions,
                                fields=[("enzo", fields[0]), 
                                        ("enzo", fields[1]),
                                        ("enzo", fields[2])])
        Acc = [cube[("enzo",fields[0])].d,
              cube[("enzo", fields[1])].d,
              cube[("enzo", fields[2])].d]

        Dil = getRotFreeField(Acc)

        Rot = np.array([Acc[0] - Dil[0],
                    Acc[1] - Dil[1],
                    Acc[2] - Dil[2]])

# safety check
        if 0:
            if (np.abs(curl2(Dil)).max() > 1e-10):
                print "whoopsie, sth went wrong in calculation the dilatational part"
                
            if (np.abs(div2(Rot)).max() > 1e-10):
                print "whoopsie, sth went wrong in calculation the rotational part"            

        RotSqr = 0.
        DilSqr = 0.

        for i in range(3):
            RotSqr += np.sum(np.abs(fftn(Rot[i]))**2.)
            DilSqr += np.sum(np.abs(fftn(Dil[i]))**2.)   


        rat= (RotSqr/(DilSqr + RotSqr))



        dpower = 0
        rpower = 0
        vpower = 0
        for dim in range(3):
            vdilhat = fftn(Dil[dim])
            vrothat = fftn(Rot[dim])
            vhat    = fftn(Acc[dim])
            dpower = vdilhat*np.conj(vdilhat) + dpower
            rpower = vrothat*np.conj(vrothat) + rpower
            vpower = vhat*np.conj(vhat) + vpower
        ffd = Filter.FourierFilter(dpower)
        ffr = Filter.FourierFilter(rpower)
        ffv = Filter.FourierFilter(rpower)
        power_1d = np.array([dpower[ffd.get_shell(bin)].sum() for bin in range(ffd.nx)])
        power_1r = np.array([rpower[ffr.get_shell(bin)].sum() for bin in range(ffr.nx)])
        power_1v = np.array([vpower[ffv.get_shell(bin)].sum() for bin in range(ffv.nx)])
        kspace=ffd.get_shell_k()
        mk = kspace[ kspace != 0].min()
        fPickle.dump(power_1d,'p42_power_1d_%s_%04d.pickle'%(car.name,frame))
        fPickle.dump(power_1r,'p42_power_1r_%s_%04d.pickle'%(car.name,frame))
        fPickle.dump(power_1v,'p42_power_1v_%s_%04d.pickle'%(car.name,frame))
        fPickle.dump(kspace,'p42_kspace_%s_%04d.pickle'%(car.name,frame))

        plt.clf()
        kmk = (kspace/mk)[1:]
        plt.plot(kmk, power_1v[1:]*kmk**(5./3),label='v')
        plt.xscale('log'); plt.xlabel(r'$k/k_{\rm{min}}$')
        plt.yscale('log'); plt.ylabel(r'$P(k) k^{5/3}$')
        outname = 'p42_power_comp_%s_%04d.pdf'%(car.name, frame)
        plt.savefig(outname)
        print outname
        plt.clf()
        plt.plot((kspace/mk)[1:], power_1d[1:],label='comp')
        plt.plot((kspace/mk)[1:], power_1r[1:],label='sol')
        plt.xscale('log') ; plt.xlabel(r'$k/k_{\rm{min}}$')
        plt.yscale('log'); plt.ylabel(r'$P_{\rm{sol}}(k), P_{\rm{comp}}(k)$')
        outname = 'p42_helm_power_%s_%04d.pdf'%(car.name,frame)
        plt.legend(loc=0)
        plt.savefig(outname)
        plt.clf()
        plt.plot((kspace/mk)[1:], power_1r[1:]/(power_1r[1:]+power_1d[1:]),label='sol/(sol+comp)')
        plt.xscale('log'); plt.xlabel(r'$k/k_{\rm{min}}$')
        plt.ylabel(r'$\chi (1-\chi?)$')
        outname = 'p42_power_ratio_%s_%04d.pdf'%(car.name,frame)
        plt.legend(loc=0)
        plt.savefig(outname)
        print outname
        plt.clf()
        plt.plot((kspace/mk)[1:], power_1v[1:],label='v')
        plt.xscale('log'); plt.xlabel(r'$k/k_{\rm{min}}$')
        plt.yscale('log'); plt.ylabel(r'$P(k)$')
        outname = 'p42_power_%s_%04d.pdf'%(car.name, frame)
        plt.savefig(outname)
        print outname

