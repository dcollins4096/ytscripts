<<<<<<< working copy: e6985a88e431 - dcollins4096: temp still not working
my_ts = tsA01 #tsy701
my_stuff = stuff
directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/rA01_rb96_110_f-/%s'
#directory = '/Users/dcollins/scratch/Paper49b_play/Eigen/y701_rb96_fft_f-_play/%s'
off_disk = {}
names={}
n2={}
n2['hx']='Bx'
n2['hy']='By'
n2['hz']='Bz'
n2['d']='density'
n2['vx']='x-velocity'
n2['vy']='y-velocity'
n2['vz']='z-velocity'
names['hx']='Bx_16.h5'
names['hy']='By_16.h5'
names['hz']='Bz_16.h5'
names['d']='density_16.h5'
names['vx']='x-velocity_16.h5'
names['vy']='y-velocity_16.h5'
names['vz']='z-velocity_16.h5'
#['d','density'],['vx','x-velocity'],['hx','Bx'],['hy','By'],['hz','Bz'],
#                ['vz','z-velocity'],['vy','y-velocity'], ['p','GasPressure']
for field in ['hx']: #['d','vx','vy','vz','hx','hy','hz']:# field_list:
    #print( field, np.abs(my_ts.temp_means[field])-my_stuff['means'][field])

    lc = my_stuff['cubes'][field]
    oc = my_ts.temp_cubes[field]
    occ = my_ts.cubes[n2[field]]

    df = h5py.File(directory%names[field])
    dc_all = df[names[field]][:]
    dcc = dc_all.swapaxes(0,2)
    dc = dc_all.swapaxes(0,2)[:16,:16,:16]
    df.close()

    ff = h5py.File(directory%('DD0000/data0000.cpu0000'))
    fc_all = ff['Grid00000001']['BxF'][:]
    fcc = fc_all.swapaxes(0,2)
    fc = fc_all.swapaxes(0,2)[:16,:16,:16]
    ff.close()
    gf = h5py.File(directory%('DD0000/data0000.cpu0000'))
    gc_all = gf['Grid00000001']['Bx'][:]
    gcc = gc_all.swapaxes(0,2)
    gc = gc_all.swapaxes(0,2)[:16,:16,:16]
    ff.close()
    
    bar = 0.# np.mean(oc)
    print( field, dif)
    print("lc  ",lc[:5,8,8]-bar)
    print("occ ",occ[:5,8,8]-bar)
    print("dcc ",dcc[:5,8,8]-bar)
    print( "dc - oc %0.2e"%(np.abs(dc-oc).sum()))
    print( "dc - lc %0.2e"%(np.abs(dc-lc).sum()))
    print( "dcc - occ %0.2e"%(np.abs(dcc-occ).sum()))
||||||| base

if 'ef' not in dir():
    execfile('go')
reload(taxi)

sys.path.append('/Users/dcollins/local-other-2018-01-05/cmbtools_nofits')
import turb_quan
reload(turb_quan)
import p49_fields
reload(p49_fields)
car=taxi.taxi('r22')
if 1:
    car.derived_fields['QU'] = p49_fields.add_QU
if 0:
    car.load()
    ad=car.ds.all_data()
    stat(ad['Ux'])
    car.plot()
if 1:
    qb = turb_quan.quan_box(car)
    qb.make_frbs(0, ['z'])
    qb.QUEB(0)
if 1:
    qb = turb_quan.quan_box(car)
    qb.GetQUEB(0)
if 1:
    plt.clf()
    for field in 'QUEB':
        for f in qb.QUEBarr[field]:
            oot = "%s_%s"%(car.name,f.split('/')[-1])
            plt.clf()
            plt.imshow(qb.QUEBarr[field][f],interpolation='nearest',origin='lower')
            plt.title(oot)
            plt.savefig(oot+".png")
            print(oot)
    

if 0:
    """
    Jordan Mirocha
    cosmorec chluba & thomas
    Wouthuysen Field Effect (Vowt Howsen) Prichartd & Furlanetto 2006
    Furlaneto 2006
    Mesinger+2011
    deOlivera-Costa + 2008 sky model
    loco.lab.asu.edu/download

"""
if 0:
    axis='z'
    field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
    field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]

    def _Q_local(field,data):
        """This function calculates the Stokes Parameter "Q" along an axis x, y, or z.
        Makes use of the depolarization factor "epsilon" using a power exponent.
        """
        n = data['density']
        B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0

        #epsilon = np.ones(data['density'].shape)
        #epsilon[ n <= n0 ] = (n.v)[ n <= n0 ]  
        #epsilon[ n > n0 ]  = (n0**(1-p) * n.v**p)[ n > n0 ]   
        epsilon = n
        out = ( epsilon * (((data[field_horizontal])**2.0) - ((data[field_vertical])**2.0))/B_sq )
        out[B_sq == 0] = 0
        return out
    car.load()
    ad=car.ds.all_data()
    b=_Q_local(None,ad)
    stat(b)
=======

if 'ef' not in dir():
    execfile('go')
reload(taxi)
import turb_quan
reload(turb_quan)

import xtra_energy_fields
reload(xtra_energy_fields)
car = taxi.taxi('ab25')
car.load()
field = 'eng_x'
#car.ds.add_field(field,**xtra_energy_fields.field_args[field])
xtra_energy_fields.dave_add_field(car.ds) #, field_name = field)
print(car.ds.index.grids[0][field])
#qb = turb_quan.quan_box(car)
#qb.load()
field_args={}
>>>>>>> merge rev:    d6fe13cfeff2 - dcollins4096: tweaks
