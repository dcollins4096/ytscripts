import dsdiff_helpers
import dsdiff
reload(dsdiff)

#basedir = '/Users/dcollins/scratch/EnzoProjects/E21_AdiabaticCorners'
basedir = '/Users/dcollins/scratch/EnzoProjects/E12/E12a_AdiabaticExpansion'
#basedir = '/Users/dcollins/scratch/EnzoProjects/E12/E12c_Zeldovich'
e01 = 'e01_ae_ppm'
e02 = 'e02_ae_ded'
e03 = 'e03_ae_ct'
e04 = 'e04_ae_ded_old'
e06 = 'e06_ae_amr_ded'
e08 = 'e08_ded_v'
e08b = 'e08b_repeat'
e09 = 'e09_ct_v'
e10 = 'e10_PR_ded'
e10b= 'e10b_repeat'
e10c= 'e10c_like_e11'
e10d= 'e10d_repeat_again'
e11 = 'e11_PR_ct'
e12 = 'e12_PR_ded_amr'
e14 = 'e14_ded_PostEnzo15'
e20 = 'e20_ded_final'
e21 = 'e21_ct_final'
e22 = 'e22_ded_orig'
sim1 = e10
#sim_list = [e06, e12]
#sim_list = [e02,e03] #,e08b]
sim_list = [e10d, e20,e21,  e22]

#a06 = 'a06_ded_hack'
#a08 = 'a08_ded_no_hack'
#a11 = 'a11_ct_goofoof'
#sim_list = [a06,a08,a11]
#MHDSimList = [e02,e04, e03, e08, e09,e11,e10]
NonMHDSimList = [e01]
color_dict = {e02:'k',e04:'c', e03:'m', e08:'c', e09:'m',e11:'g',e10:'r',e12:'r',e10b:'r'}
color_dict = {}
def random_color():
    return [random.random(),random.random(),random.random()]
#sim1 = 'a08_post_sbc'
if 'frame' not in dir():
    frame = 3

def prnt(normlist):
    keys = sorted(normlist.keys())
    for ke in keys:
        values = normlist[ke]
        print "%15s"%ke, "%9.2e "*len(values)%tuple(values)

#dir_name1 = '%s/%s'%(basedir,sim1)
#dir_name2 = '%s/%s'%(basedir,sim2)
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
        E_list_plot = {}
        E_norm = {}
        B_norm = {}
        B_list={}
        readit = True
    plt.clf()
    for sim in sim_list:
        if readit:
            a_list[sim] = []
            E_list[sim] = []
            B_list[sim] = []
            for frame in [0,2,4,10, 12,16]: #range(17):

                if 0:
                    dir_name = '%s/%s'%(basedir,sim)
                    ds_name = dsdiff_helpers.get_ds_name(dir_name,frame)
                    ds = yt.load(ds_name)
                    this_z = ds['CosmologyCurrentRedshift']
                    this_a = 1./(1.+this_z)
                    ad = ds.all_data()
                    this_e = ad['TotalEnergy']
                    this_e_avg = (this_e*ad['cell_volume']).sum()/ad['cell_volume'].sum()
                    uniformity_error = (np.abs(this_e - this_e[0])).max() 
                    if uniformity_error > 1e-9:
                        print "WARNING non uniform", uniformity_error
                    if sim not in NonMHDSimList: #'e02_ae_ded':
                        Bx = ad['Bx'].in_units('code_magnetic')
                        By = ad['By'].in_units('code_magnetic')
                        Bz = ad['Bz'].in_units('code_magnetic')
                        this_b = 0.5*(Bx*Bx+By*By+Bz*Bz)
                        if np.abs(this_b-this_b[0]).sum() > 1e-9:
                            print "WARNING: nonuniform B"
                        B_list[sim].append(this_b[0])

                if 1:
                    dir_name = '%s/%s'%(basedir,sim)
                    grid_name = dsdiff_helpers.find_grid_filename(dir_name,frame,1)
                    fptr = h5py.File(grid_name)
                    grid = 1
                    dir_name = '%s/%s'%(basedir,sim)
                    ds_name = dsdiff_helpers.get_ds_name(dir_name,frame)
                    ds = yt.load(ds_name)
                    this_z = ds['CosmologyCurrentRedshift']
                    this_a = 1./(1.+this_z)
                    this_e = fptr['Grid%08d'%grid]['TotalEnergy'][:].flatten()
                    this_e_avg = this_e.sum()/this_e.size
                    uniformity_error = (np.abs(this_e - this_e[0])).max() 
                    if uniformity_error> 1e-9:
                        print "WARNING non uniform", uniformity_error
                    if sim not in NonMHDSimList: #'e02_ae_ded':
                        Bx = fptr['Grid%08d'%grid]['Bx'][:].flatten()
                        By = fptr['Grid%08d'%grid]['By'][:].flatten()
                        Bz = fptr['Grid%08d'%grid]['Bx'][:].flatten()
                        this_b = 0.5*(Bx*Bx+By*By+Bz*Bz)
                        if np.abs(this_b-this_b[0]).sum() > 1e-9:
                            print "WARNING: nonuniform B"
                        B_list[sim].append(this_b[0])

                a_list[sim].append(this_a)
                E_list[sim].append(this_e_avg)
        a_list[sim] = nar(a_list[sim])
        B_list[sim] = nar(B_list[sim])
        E_list[sim] = nar(E_list[sim])
        a0=a_list[sim][0]
        if 1:
            E_list_plot[sim] = E_list[sim]-B_list[sim]
        else:
            E_list_plot[sim] = copy.copy(E_list[sim])
        E0=E_list_plot[sim][0]

        E_norm[sim]=(E_list_plot[sim]/E_list_plot[sim][0])**-0.5/a_list[sim]*a_list[sim][0]-1.0
        #E_norm[sim] += {e01:0,e02:1e-2, e03:2e-2, e04:3e-2}[sim]
        color =  color_dict.get(sim,random_color())
        plt.plot(a_list[sim],E_norm[sim], c=color,marker='+',label=r'$E^{-1/2}/a\ %s$'%(sim[0:3]),
                 linestyle = "--")
        #plt.plot(a_list[sim],a_list[sim], c=ecolor[sim],marker='+',label=r'$E^{-1/2}/a\ %s$'%(sim[0:3]))
        if sim not in NonMHDSimList:
            B_norm[sim] = (B_list[sim]/B_list[sim][0])**-1./a_list[sim]*a_list[sim][0]-1.0
            plt.plot(a_list[sim],B_norm[sim], c=color, label=r'$B^{-1}/a\ %s$'%(sim[0:3]))
        #plt.plot( a_list[sim], B_list[sim],c={e02:'k',e04:'c',e03:'m'}[sim],label='B %s'%(sim[0:3]))
        print "E/B (%0.2e/%0.2e) %0.2e"%(E_list_plot[sim][0], B_list[sim][0], E_list_plot[sim][0]/B_list[sim][0])
        
    n_sims = len(sim_list)
    sim_string = "_%s"*n_sims%tuple(sim_list)
    outname = 'AdiatabicExpansionPlot%s.pdf'%(sim_string) #'%s_ae.pdf'%(sim)
    outname = 'Adiabatic.pdf'
    #plt.ylim(-1,1)
    #plt.yscale('log')
    plt.ylabel('ShouldBeZero')
    plt.xlabel(r'$a$')
    plt.legend(loc=2)
    plt.savefig(outname)
    print "a list"
    prnt(a_list)
    print "E norm"
    prnt(E_norm)
    print "B norm"
    prnt(B_norm)
    print outname
            

