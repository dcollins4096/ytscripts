
plt.clf()
delta = T2-T4
wtfmate=np.abs(delta[:,:,16])
p=plt.imshow( np.log10(wtfmate), interpolation='nearest',origin='lower')
plt.colorbar(p)
outname = "p35_dt_%s_%04d_%s"%(c05x.name, c05x.frames[0], 'delta_velocity')
plt.savefig(outname)
print outname
