from yt import YTArray as yta
import yt.units as units
G = units.gravitational_constant
c = units.speed_of_light

class snapshot():
    def __init__(self, npoints=0):
        self.T = -1
        self.dT = -1
        self.WhateverThisIs = -1
        self.data = np.zeros([npoints,78])-1.234
        self.slice = slice(None)
    def __getitem__(self,item_in):
        #stuff map (n_column, units)
        c_columns=[  'n', 'p', 'd', 'he3', 'he4', 'li7', 'be7', 'be8', 'b8', 'c12', 'c13', 'n13', 'n14', 'n15', 'o15', 'o16', 'o17',
          'o18', 'f17', 'f18', 'f19', 'ne20', 'ne21', 'ne22', 'na22', 'na23', 'na24', 'mg24', 'mg25', 'mg26', 'al25',
          'al26', 'al27']
        k_columns=['id', 'M(r)', 'dm', 'r[cm]', 'v[cm/s]', 'T[K]', 'rho[g/ccm]', 'U[erg]', 'dU/dTemp', 'p[erg]', 'dp/dTemp', 'dEnuc/dt']

        units={'M(r)':'Msun','r[cm]':'cm','T[K]':'K','rho[g/ccm]':'g/cm**3', 'dm':'Msun'}
        alias={'M':'M(r)','r':'r[cm]','rho':'rho[g/ccm]'}
        
        if item_in in alias:
            item = alias[item_in]
        else:
            item = item_in
        if item in k_columns:
            nItem = k_columns.index(item)
        if item in c_columns:
            nItem = 2*c_columns.index(item) + len(k_columns)
        these_units = units.get(item,None)
        output = yta(self.data[:,nItem],these_units)
        return output[self.slice]

    def plot(self,a,b,**kwargs):
        plt.plot( self[a], self[b],**kwargs )
        plt.xlabel(a)
        plt.ylabel(b)
        outname = "%s_%s.png"%(a,b)
        print(outname)
        return outname 
    def error_check(self):
        n_data = ( np.abs( self.data + 1.234 ) < 1e-6 ).sum()
        print("Number not filled: %d"%( n_kin))

def consume(fname, snapshots=[], nstop=-1,npoints = 910, nskip=None):
    fptr = open(fname,'r')
    timestamp_line = False
    label_line = False
    data_line=True
    #might not work for all outputs.
    this_snap = snapshot(npoints=npoints)
    snapshots.append(this_snap)
    for nline,line in enumerate(fptr):
        if nskip is not None:
            if nline < nskip:
                continue
        spl = no_whites(line.split(' '))
        if len(spl) not in [1, 3, 13, 14, 15, 45, 78]:
            print("Parse error: line %d has %d elements"%(nline,len(spl)))
            print(line)
            print(" (you really should make this an exception...)")
            snapshots.append(line)
            return snapshots
        if len(spl) == 1:
            continue
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
        if data_line:
            v = nar(spl).astype('float')
            nelement=int(v[0])-1
            this_snap.data[nelement,:] = v

        if nelement == npoints-1:
            data_line=False
            head_line=True


        if nline == nstop:
            return snapshots
    return snapshots

fname = '/home/dcollins/RESEARCH2/Paper52_NewWD/2018-10-30-moreIC/wd.struc_at_runaway'

snaps=consume(fname, nskip = 4 )
ts = snaps[0]
#print(ts['M'])
def mid(arr):
    return 0.5*(arr[1:]+arr[:-1])
def delta(arr):
    return (arr[1:]-arr[:-1])
dr = np.zeros_like(ts['r'])
dr[:-1] = delta(ts['r'])
dm = (ts['rho']*4*np.pi*ts['r']**2*dr).in_units('msun')

