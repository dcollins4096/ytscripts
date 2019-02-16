import get_all_quantities as gaq
reload(gaq)

if 'car' not in dir():
    car = taxi.load('p49d_za01')
    car.load(0)
    quan = gaq.all_quan_from_taxi(car)
#plt.clf()
#plt.plot(quan['time'],quan['density_std'],marker="*")
#plt.savefig('../PigPen/%s_time_density_std.png'%car.outname)
v2 = np.sqrt(quan['vx_std']**2+quan['vy_std']**2+quan['vz_std']**2)
b2 = np.sqrt(quan['bx_std']**2+quan['by_std']**2+quan['bz_std']**2+
             quan['bx_avg']**2+quan['by_avg']**2+quan['bz_avg']**2)/np.sqrt(np.pi*4)

plt.clf()
plt.plot(quan['time'],v2/b2,marker="*")
plt.savefig('../PigPen/%s_time_alfmach.png'%car.outname)
plt.clf()
plt.plot(quan['time'],v2,marker="*")
plt.savefig('../PigPen/%s_time_mach.png'%car.outname)
plt.clf()
plt.plot(quan['time'],b2,marker="*")
plt.savefig('../PigPen/%s_time_b.png'%car.outname)

