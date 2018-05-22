
if 'ef' not in dir():
    execfile('go')
reload(taxi)

sys.path.append('/Users/dcollins/local-other-2018-01-05/cmbtools_nofits')
import turb_quan
reload(turb_quan)
import p49_fields
reload(p49_fields)
car=taxi.taxi('r22')
car.outname='r22_cmbtools_0.1_'
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
            oot = "%s_%s"%(car.outname,f.split('/')[-1])
            plt.clf()
            f = plt.imshow(qb.QUEBarr[field][f],interpolation='nearest',origin='lower')
            plt.colorbar(f)
            print("did a color")
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
