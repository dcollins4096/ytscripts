
from yt.analysis_modules.ppv_cube.api import PPVCube
import yt.units as u
def fake_temp(field,data):
    return data.ds.arr(np.ones_like(data['density'].v)*10,'K')

yt.add_field('temperature',function=fake_temp,units='K')

if 0:

    simname = 'B02'
    resolution = 256
    frame = 30
    dirname = '/scratch1/dcollins/Paper08/%s/%d/RS%04d/restart%04d'%(simname,resolution,frame,frame)

    ds = yt.load(dirname)
    ds.index.grids[0]['temperature']
    L = [0,0,1]
    cube = PPVCube(ds, L, "density", (-150.,150.,50,"code_velocity"), dims=200, method="sum")

if 0:
    for n,v in enumerate(cube.vmid):
#       if n not in [20,21]:
#           continue
        plt.clf()
        plt.imshow(np.log10(cube[:,:,n].v),cmap='jet',origin='lower',interpolation='nearest')
        plt.text(0,0,r'$%0.2f \rm{km}/\rm{s}$'%v)
        cb=plt.colorbar()
        cb.cmap.set_under('w')
        cb.cmap.set_bad('w')
        outname = 'p14b_tmp_slice%04d.png'%n
        plt.savefig(outname)
        print outname

if 1:
    plt.clf()
    cube_read,cube_ds,vel = at.get_ppv('x',prefix='ppv_1', ppv_args={'dims':256,'velocity_bounds':(-35,35,14,'code_velocity')})
    plt.imshow(np.log10(np.sum(cube_read[:,:,:],axis=2)),cmap='jet',origin='lower',interpolation='nearest')
    cb=plt.colorbar()
    cb.cmap.set_under('w')
    cb.cmap.set_bad('w')
    outname = 'p14b_tmp_integral'
    plt.savefig(outname)
    print outname

