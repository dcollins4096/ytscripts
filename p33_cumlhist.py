
if 'fltc' not in dir():
    fltc=taxi.fleet(['aj25', 'aj26']) #,'aj22'])

    def _dbg(field,data):
        return data['DebugField']
    yt.add_field('SFdbg',function=_dbg,take_log=False)

    ef('p33_sfhunt.py')

    fltc['frames']=[18]
    fltc['fields'] = ['mjeans','bmass']
    fltc('car.outname = car.name+"_t2"')
    fltc['callbacks']=[]
fltc('car.fill(18)')
g1 = fltc[0].ds.index.grids[-1]
g2 = fltc[1].ds.index.grids[-1]
plt.clf()
b1 = g1['bmass'].flatten().v
b2 = g2['bmass'].flatten().v
j1 = g1['mjeans'].flatten().v
j2 = g2['mjeans'].flatten().v
frm1 = g1['sf_jeans']
frm2 = g2['sf_jeans']

plt.clf()
plt.plot( sorted(b1), np.cumsum(sorted(b1)), label='b h0', c='r')
plt.plot( sorted(b2), np.cumsum(sorted(b2)), label='b h6', c='g')
plt.plot( sorted(j1), np.cumsum(sorted(j1)), label='j h0', c='r', linestyle='--')
plt.plot( sorted(j2), np.cumsum(sorted(j2)), label='j h6', c='g', linestyle='--')
plt.xlabel('M')
plt.ylabel(r'$\int_0^M M dM')
plt.legend(loc=0)
plt.xscale('log'); plt.yscale('log')
plt.savefig('bmass_%s.pdf'%fltc.allnames())

ok1 = b1 > j1
ok2 = b2 > j2


stat(b1,'b1')
stat(b2,'b2')
stat(j1,'j1')
stat(j2,'j2')
