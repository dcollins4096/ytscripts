import get_all_quantities as gaq
reload(gaq)


flt = taxi.fleet(['za01','zb01','zc01_quan','zd01_quan','ze01_quan'])
def plot_quan(car):
    if type(car) is str:
        car = taxi.load(car)
    
    quan = gaq.all_quan_from_taxi(car)
    v2 = np.sqrt(quan['vx_std']**2+quan['vy_std']**2+quan['vz_std']**2)
    b2 = np.sqrt(quan['bx_std']**2+quan['by_std']**2+quan['bz_std']**2+
                 quan['bx_avg']**2+quan['by_avg']**2+quan['bz_avg']**2)/np.sqrt(np.pi*4)

    plt.clf()
    plt.plot(quan['time'],np.sqrt(v2/b2),marker="*")
    plt.savefig('../PigPen/%s_time_alfmach.png'%car.outname)
    plt.xlabel('time')
    plt.ylabel('Alfven Mach $v_{\rm{rm}}/(b_{\rm{rms}}/\bar{rho}')
    plt.clf()
    plt.plot(quan['time'],v2,marker="*")
    plt.savefig('../PigPen/%s_time_mach.png'%car.outname)
    plt.clf()
    plt.plot(quan['time'],b2,marker="*")
    plt.savefig('../PigPen/%s_time_b.png'%car.outname)
    return {'quan':quan,'car':car}

