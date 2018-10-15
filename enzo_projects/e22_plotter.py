
execfile('go_lite')
import yt
import dsdiff_helpers
base_dir = '/scratch1/dcollins/EnzoProjects/E22_DEF_NonCosmo/'
plot_sequence_name = 'E22_Turb_full'
thelist = ['f01_ppm_nodef', 'f02_ct_nodef', 'f03_ded_nodef']
thelist = ['f01c_ppm_probM1']
thelist = ['f01_ppm_nodef',  'f03_ded_nodef', 'f05_ct_newversion_newparam','f01b_ppm_problem0']
thelist += [ 'f06_ppm_hll']
plot_sequence_name = 'E22_Turb_g'
field = ("enzo","TotalEnergy"); ofield = 'enzo_TotalEnergy'
field_list = ['temperature','density']
for frame in [1]+range(10,60,10):# range(50): #[0,1,2,5,10,15,20]: #framelist.get(dirname, default_frames):
    plt.clf()
    for dirname in thelist:
        print "frame", frame
        short_name = dirname[0:6]
        setname = dsdiff_helpers.get_ds_name("%s/%s"%(base_dir,dirname),frame)
        print setname
        ds = yt.load(setname)
        for field in field_list:
            print yt.ProjectionPlot(ds,0,field).save('%s_%s_%04d'%(plot_sequence_name, short_name, frame))


