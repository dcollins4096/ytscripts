
aq42 = taxi.taxi(directory='/scratch/00369/tg456484/Paper42_NewAK/aq42_m9_drive0.5_512_p59_ppm',name='aq42',frames=range(90),fields=['density'])
car = aq42
frame = 90
frames = [40,90]
power_1d = []
power_1r = []
power_1v = []
kspace = []
for frame in frames:
    power_1d.append(fPickle.load('p42_power_1d_%s_%04d.pickle'%(car.name,frame)))
    power_1r.append(fPickle.load('p42_power_1r_%s_%04d.pickle'%(car.name,frame)))
    power_1v.append(fPickle.load('p42_power_1v_%s_%04d.pickle'%(car.name,frame)))
    kspace.append(fPickle.load('p42_kspace_%s_%04d.pickle'%(car.name,frame)))
mk = kspace[0][ kspace != 0].min()

plt.clf()
lines=['--',':']
frame_str = "_%04d"*len(frames)%tuple(frames)
for n in range(len(frames)):
    kmk = (kspace[n]/mk)[1:]
    plt.plot(kmk, power_1v[n][1:]*kmk**(5./3),label='v t=%d'%frames[n],linestyle=lines[n])
plt.xscale('log'); plt.xlabel(r'$k/k_{\rm{min}}$')
plt.yscale('log'); plt.ylabel(r'$P(k) k^{5/3}$')
outname = 'p42_power_comp_%s%s.pdf'%(car.name, frame_str)
plt.legend(loc=0)
plt.savefig(outname)
print outname
plt.clf()
for n in range(len(frames)):
    plt.plot((kspace[n]/mk)[1:], power_1d[n][1:],label='comp_%04d'%frames[n], linestyle=lines[n])
    plt.plot((kspace[n]/mk)[1:], power_1r[n][1:],label='sol_%04d'%frames[n], linestyle=lines[n])
plt.xscale('log') ; plt.xlabel(r'$k/k_{\rm{min}}$')
plt.yscale('log'); plt.ylabel(r'$P_{\rm{sol}}(k), P_{\rm{comp}}(k)$')
outname = 'p42_helm_power_%s_%s.pdf'%(car.name,frame_str)
plt.legend(loc=0)
plt.savefig(outname)
print outname
plt.clf()
for n in range(len(frames)):
    plt.plot((kspace[n]/mk)[1:], power_1r[n][1:]/(power_1r[n][1:]+power_1d[n][1:]),label='sol/(sol+comp) %04d'%frames[n], linestyle=lines[n] )
plt.xscale('log'); plt.xlabel(r'$k/k_{\rm{min}}$')
plt.ylabel(r'$\chi (1-\chi?)$')
outname = 'p42_power_ratio_%s%s.pdf'%(car.name,frame_str)
plt.legend(loc=0)
plt.savefig(outname)
print outname
plt.clf()
for n in range(len(frames)):
    plt.plot((kspace/mk)[1:], power_1v[1:],label='v %04d'%frames[n], linestyle=lines[n])
plt.xscale('log'); plt.xlabel(r'$k/k_{\rm{min}}$')
plt.yscale('log'); plt.ylabel(r'$P(k)$')
outname = 'p42_power_%s%s.pdf'%(car.name, frame_str)
plt.savefig(outname)
print outname

