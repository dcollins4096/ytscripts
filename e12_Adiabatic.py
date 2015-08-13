import dsdiff_helpers
import dsdiff
reload(dsdiff)

#basedir = '/Users/dcollins/scratch/EnzoProjects/E21_AdiabaticCorners'
basedir = '/Users/dcollins/scratch/EnzoProjects/E12/'
e01 = 'e01_ae_ppm'
e02 = 'e02_ae_ded'
e03 = ''
e04 = 'e04_ae_ded_old'

sim1 = e04
sim2 = e02    
MHDSimList = [e02,e04]
#sim1 = 'a08_post_sbc'
if 'frame' not in dir():
    frame = 3

dir_name1 = '%s/%s'%(basedir,sim1)
dir_name2 = '%s/%s'%(basedir,sim2)
if 0:
    ds_name = dsdiff_helpers.get_ds_name(dir_name,frame)
    print ds_name
    ds = yt.load(ds_name)
    ad = ds.all_data()
    den = ad['density']
    stat(den.in_units('code_density').v-1.0,'den')
    proj = yt.ProjectionPlot(ds,0,'density')
    proj.annotate_grids()
    outname = sim1+'%04d'%frame
    print proj.save(outname)

if 0:
    dir1 = dir_name1
    dir2 = dir_name2
    fields = adiabatic_hydro
    frames=[frame]
    grids=[1]
    ud = dsdiff.udiff(dir1,dir2,frames=frames,grids=grids, fields=['TotalEnergy'])
    ud()
    #ud(shift=[1.0,0.0,0.0])

if 1:
    if 'a_list' not in dir():
        a_list = {}
        E_list = {}
        E_norm = {}
        B_norm = {}
        B_list={}
    plt.clf()
    for sim in [sim1,sim2]:
        if readit:
            a_list[sim] = []
            E_list[sim] = []
            B_list[sim] = []
            for frame in [0,2,4,10, 12,16]: #range(17):
                dir_name = '%s/%s'%(basedir,sim)
                ds_name = dsdiff_helpers.get_ds_name(dir_name,frame)
                ds = yt.load(ds_name)
                this_z = ds['CosmologyCurrentRedshift']
                this_a = 1./(1.+this_z)
                ad = ds.all_data()
                this_e = ad['TotalEnergy']
                if (np.abs(this_e - this_e[0])).max() > 1e-9:
                    print "WARNING non uniform"
                if sim in MHDSimList: #'e02_ae_ded':
                    Bx = ad['Bx'].in_units('code_magnetic')
                    By = ad['By'].in_units('code_magnetic')
                    Bz = ad['Bz'].in_units('code_magnetic')
                    this_b = 0.5*(Bx*Bx+By*By+Bz*Bz)
                    if np.abs(this_b-this_b[0]).sum() > 1e-9:
                        print "WARNING: nonuniform B"
                    B_list[sim].append(this_b[0])
                a_list[sim].append(this_a)
                E_list[sim].append(this_e[0])
        a_list[sim] = nar(a_list[sim])
        B_list[sim] = nar(B_list[sim])
        E_list[sim] = nar(E_list[sim])
        a0=a_list[sim][0]
        E0=E_list[sim][0]

        E_norm[sim]=(E_list[sim]/E_list[sim][0])**-0.5/a_list[sim]*a_list[sim][0]-1.0
        E_norm[sim] += {e01:0,e02:1e-2, e03:2e-2, e04:3e-2}[sim]
        ecolor={sim1:'r',sim2:'g'}
        plt.plot(a_list[sim],E_norm[sim],
                c=ecolor[sim],marker='+')
        if sim in MHDSimList:
            B_norm[sim] = (B_list[sim]/B_list[sim][0])**-1./a_list[sim]*a_list[sim][0]-1.0
            plt.plot(a_list[sim],B_norm[sim],
                     c={e02:'k',e04:'c'}[sim])
    outname = 'thing.pdf' #'%s_ae.pdf'%(sim)
    plt.savefig(outname)
    print outname
            

