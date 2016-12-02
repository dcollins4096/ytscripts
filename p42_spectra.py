import p42_pg


aq32 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq32_m2.9_drive1_32',name='aq32',frames=range(0,42,2),fields=['density'])
aq31 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq31_m2.9_drive0.5_32',name='aq31',frames=range(0,42,2),fields=['density'])
aq16 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq16_m2.9_drive0_noamr_128',name='aq15',frames=range(17),fields=['density'])
car = aq31
frame = 40
car.fill(frame)
ds=car.ds

fields = ['%s-acceleration'%s for s in 'xyz']
#fields = ['DrivingField%s'%s for s in '123']
cube = ds.covering_grid(0, left_edge=ds.domain_left_edge,
                        dims=ds.domain_dimensions,
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



import fourier_tools.fourier_filter as Filter
dpower = 0
rpower = 0
for dim in range(3):
    vdilhat = fftn(Dil[dim])
    vrothat = fftn(Rot[dim])
    dpower = vdilhat*np.conj(vdilhat) + dpower
    rpower = vrothat*np.conj(vrothat) + rpower
ffd = Filter.FourierFilter(dpower)
ffr = Filter.FourierFilter(rpower)
power_1d = np.array([dpower[ffd.get_shell(bin)].sum() for bin in range(ffd.nx)])
power_1r = np.array([rpower[ffr.get_shell(bin)].sum() for bin in range(ffr.nx)])
kspace=ffd.get_shell_k()
mk = kspace[ kspace != 0].min()
plt.clf()
plt.plot((kspace/mk)[1:], power_1d[1:],label='dil')
plt.plot((kspace/mk)[1:], power_1r[1:],label='rot')
plt.xscale('log')
plt.yscale('log')
outname = 'p42_power_%s_%04d.pdf'%(car.name,frame)
plt.legend(loc=0)
plt.savefig(outname)
plt.clf()
plt.plot((kspace/mk)[1:], power_1r[1:]/(power_1r[1:]+power_1d[1:]),label='rot/(rot+dil)')
plt.xscale('log')
outname = 'p42_power_ratio_%s_%04d.pdf'%(car.name,frame)
plt.legend(loc=0)
plt.savefig(outname)
print outname



