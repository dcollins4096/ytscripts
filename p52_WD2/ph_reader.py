from yt import YTArray as yta
import yt.units as units
G = units.gravitational_constant
c = units.speed_of_light

class snapshot():
    def __init__(self, npoints=0):
        self.T = -1
        self.dT = -1
        self.WhateverThisIs = -1
        self.kinematics = np.zeros([npoints,14])-1.234
        self.chem = np.zeros([npoints,15])-1.234
        self.slice = slice(None)
    def __getitem__(self,item):
        #stuff map (n_column, units)
        Kmap={'id':(0,None),
                  'radius':(1,'cm'),
                  'r':(1,'cm'),
                  'm(r)':(2,'Msun'),
                  'mass':(2,'Msun'),
                  'density':(3,'g/cm**3'),
                  'velx':(4,'cm/s'),
                  'flux':(5,None),
                  'fluxc':(6,None),
                  'temp':(7,'K'),
                  'T':(7,'K'),
                  'kappa':(8,None),
                  'tau':(9,None),
                  'gammac':(10,None),
                  'gammae':(11,None),
                  'sum(Enuc)':(12,None)}
        Cmap = {'No.':(0,None),
                'he4':(1,None),
                'c12':(2,None),
                'o16':(3,None),
                'ne20':(4,None),
                'mg24':(5,None),
                'si28':(6,None),
                's32':(7,None),
                'ar36':(8,None),
                'ca40':(9,None),
                'ti44':(10,None),
                'cr48':(11,None),
                'fe52':(12,None),
                'ni56':(13,None),
                'zn60':(14,None)}
        if item in Kmap:
            output = yta(self.kinematics[:,Kmap[item][0]],Kmap[item][1])
        if item in Cmap:
            output = yta(self.chem[:,Cmap[item][0]],Cmap[item][1])
        return output[self.slice]

    def plot(self,a,b,**kwargs):
        plt.plot( self[a], self[b],**kwargs )
        plt.xlabel(a)
        plt.ylabel(b)
        outname = "%s_%s.png"%(a,b)
        print(outname)
        return outname 
    def error_check(self):
        n_chem = ( np.abs( self.chem + 1.234 ) < 1e-6 ).sum()
        n_kin = ( np.abs( self.kinematics + 1.234 ) < 1e-6 ).sum()
        print("Chem not filled: %d"%( n_chem))
        print("Kine not filled: %d"%( n_kin))
    def line_kinematic(self,line):
        vals = self.split(line," ")
        print(len(vals))

def consume(fname, snapshots=[], nstop=-1,npoints = 910):
    fptr = open(fname,'r')
    timestamp_line = True
    label_line = False
    kinematic_line=False
    chem_line=False
    for nline,line in enumerate(fptr):
        spl = no_whites(line.split(' '))
        if len(spl) not in [3, 13, 14, 15]:
            print("Parse error: line %d has %d elements"%(nline,len(spl)))
            print(line)
            print(" (you really should make this an exception...)")
            return snapshots
        if timestamp_line:
            if len(spl) != 3:
                print("line error: expected timestamp line.  Got:")
                print(line)
                return snapshots
            this_snap = snapshot(npoints=npoints)
            snapshots.append(this_snap)

            this_snap.T = float(spl[0])
            this_snap.dT = float(spl[1])
            this_snap.WhateverThisIs = float(spl[2])
            timestamp_line=False
            label_line=True
            continue
        if label_line:
            my_labels=spl
            label_line=False
            if my_labels[0] == 'id':
                kinematic_line=True
            else:
                chem_line = True
            continue
        if kinematic_line:
            v = nar(spl).astype('float')
            nelement=int(v[0])-1
            this_snap.kinematics[nelement,:] = v
        elif chem_line:
            v = nar(spl).astype('float')
            nelement=int(v[0])-1
            this_snap.chem[nelement,:] = v

        if nelement == npoints-1:
            if kinematic_line:
                kinematic_line=False
                label_line=True
            else:
                chem_line=False
                timestampe_line=True


        if nline == nstop:
            return snapshots
    return snapshots

fname = '/home/dcollins/RESEARCH2/Paper52_NewWD/2018-10-15-ICs/snia_first_snap_only_2018_10_15.plt'
snaps=consume(fname)
ts = snaps[0]
ts.slice=slice(5,None)
plt.clf()
def mid(arr):
    return 0.5*(arr[1:]+arr[:-1])
def delta(arr):
    return (arr[1:]-arr[:-1])

plt.clf()
if 1:
    ts.slice=slice(5,None)
    dPdr = -G*ts['density']*ts['r']**-2*ts['mass']
    rho = ts['density']
    rmid = mid(ts['r'])
    drhodr= delta(rho)/rmid
    ratio = (drhodr/mid(dPdr)).in_units('g/erg')
if 0:
    ts.slice=slice(5,None)
    dr=np.zeros_like(ts['r'])
    dr[0:-1] = delta( ts['r'])
    dPdr = -G*ts['density']*ts['r']**-2*ts['mass']
    rho = ts['density']
    rmid = mid(ts['r'])
    drhodr= delta(rho)/rmid
    P = np.cumsum(dP).in_units('erg/cm**3')
    plt.clf()
    plt.plot(rmid, drhodr,label=r'$d\rho/dr$')
    plt.plot(ts['r'],dPdr, label=r'$dP/dr$')
    plt.legend(loc=0)
    plt.savefig('p52_nK.png')
    plt.clf()
    ratio = (drhodr/mid(dPdr)).in_units('g/erg')

if 0:
    outname = 'p52_'+ts.plot('r','T',marker='*')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(outname)

if 0:
    ts.slice=slice(5,None)
    plt.clf()
    plt.plot(ts['r'], ts['density'])
    
    plt.savefig('p52_r_density.png')

if 0:
    """mass works like you expect."""
    dr=np.zeros_like(ts['r'])
    dr[0:-1] = delta( ts['r'])
    dm = ts['density']*ts['r']**2*dr*np.pi*4
    plt.plot(np.cumsum(dm.in_units('Msun'))-ts['mass'])
    plt.savefig('p52_dm_msun.png')

if 0:
    dr=np.zeros_like(ts['r'])
    dr[0:-1] = delta( ts['r'])
    dP = -G*ts['density']*ts['r']**-2*ts['mass']*dr
    P = np.cumsum(dP).in_units('erg/cm**3')

if 1:
    ts.slice=slice(5,None)
    plt.clf()
    #plt.plot(ts['r'],P,label='Pint')
    #P2 = ts['density']*ts['T']*units.kboltz/units.mass_hydrogen
    half = -1 #int(P.size/2)
    P += P2[half].in_units('erg/cm**3') - P[half].in_units('erg/cm**3')
    #P += np.mean(P2-P)
    plt.plot(ts['r'],P2)
    plt.plot(ts['r'],P)
    plt.yscale('log')
    plt.savefig('p52_r_P2.png')
    gamma_maybe = np.log(P)/np.log(ts['density'])
    plt.clf()
    plt.plot(ts['r'],gamma_maybe)
    plt.savefig('p52_r_gamma_maybe.png')


if 0:
    plt.plot(ts['r'],P)
    plt.xlabel(ts['r'].units)
    plt.ylabel(dP.units)
    plt.savefig('p52_r_P.png')
    

#plt.plot(ts['radius'],ts['density'])
#plt.plot(ts['radius'],ts['mass'])
#plt.plot(ts['radius'],
#plt.plot(mid(ts['r']), delta(ts['r']))
#plt.plot(ts['r'])
#plt.plot(ts['mass'])
#plt.xlabel('n')
#plt.ylabel('mass')
#plt.savefig('p52_n_mass.png')
