
if 1:
    car=aj01
    max_particle=0
    plt.clf()
    for n in particles[car.name].keys():
        p = nar(particles[car.name][n])
        L = len(p)
        plt.scatter([n]*L, p, s=0.1)
        print (p[1:]-p[:-1] - 1).sum()
    outname = 'p33_aj01_dumb_thing.png'
    plt.savefig(outname)
    print outname


if 0:
    index_list = ['aj01','aj02']
    particles={}
    particles['aj01']={}
    car = aj01
    max_id = 0
    for n in range(100):
        print "N %04d"%n,
        car.fill(n)
        if car.ds['NumberOfParticles'] > 0:
            iii = [int(mf) for mf in sorted(car.ds.all_data()['particle_index'].v)]
            particles[car.name][n]=iii
            last_p = int(iii[0])-1
            for p in iii:
                max_id = max([max_id,p])
                p = int(p)
                difference = p - last_p - 0
                for n in range(difference):
                    print "  ",
                print p,
                last_p = p
        else:
            print "--",
        print ""


