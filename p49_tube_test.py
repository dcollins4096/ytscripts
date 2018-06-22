if 'ef' not in dir():
    execfile('go')
#ef('tube.py')
import tube
reload(tube)
basedir = '/Users/dcollins/scratch/Paper49b_play/'
axis=0
coord=(0.505,0.505)
counter=0
fields=adiabatic_mhd
#fields=['density']
sims = {}
sims['a02'] =  '%s/a02_dumb_test'%basedir
sims['a03'] =  '%s/a03_transverse'%basedir
sims['r13c'] = '%s/r13c_left_fast'%basedir
sims['r14'] = '%s/r14_left_fast_test'%basedir
sims['r15'] = '%s/r15_cube_right_fast'%basedir
sims['r15b'] = '%s/r15b_test'%basedir
sims['r17'] = '%s/r17_right_python'%basedir
sims['r17b'] = '%s/r17b_right_python_test'%basedir
sims['r17c'] = '%s/r17c_check'%basedir
sims['r18'] = '%s/r18_alfven_x'%basedir
sims['r18b'] = '%s/r18b_alfven_x_fid'%basedir
sims['r19'] = '%s/r19_fast_left_rotate'%basedir
sims['r19b'] = '%s/r19b_fast_left_fid'%basedir
sims['r23a'] = '%s/r23a_fast_right_test'%basedir
sims['r23b'] = '%s/r23b_fast_left_test'%basedir
sims['r24a'] = '%s/r24a_alf_right_test'%basedir
sims['r24b'] = '%s/r24b_alf_left_test'%basedir

#fiducial square
sims['r201'] = '%s/Eigen/r201_p500_sq_f-'%basedir
sims['r202'] = '%s/Eigen/r202_p500_sq_f+'%basedir
sims['r203'] = '%s/Eigen/r203_p500_sq_s-'%basedir
sims['r204'] = '%s/Eigen/r204_p500_sq_s+'%basedir
sims['r205'] = '%s/Eigen/r205_p500_sq_a-'%basedir
sims['r206'] = '%s/Eigen/r206_p500_sq_a+'%basedir

#fiducial cos
sims['r301'] = '%s/Eigen/r301_p500_f-'%basedir

sims['r401'] = '%s/Eigen/r401_rj95_sq_f-'%basedir
sims['r402'] = '%s/Eigen/r402_rj95_sq_f+'%basedir

#fiducial cos from fft.
sims['r501'] = '%s/Eigen/r501_rj95_fft_f-'%basedir

sims['r601'] = '%s/Eigen/r601_rb96_sq_f-'%basedir
sims['r602'] = '%s/Eigen/r602_rb96_sq_f+'%basedir
sims['r603'] = '%s/Eigen/r603_rb96_sq_s-'%basedir
sims['r604'] = '%s/Eigen/r604_rb96_sq_s+'%basedir
sims['r605'] = '%s/Eigen/r605_rb96_sq_a-'%basedir
sims['r606'] = '%s/Eigen/r606_rb96_sq_a+'%basedir

sims['r701'] = '%s/Eigen/r701_rb96_fft_f-'%basedir

#rotated 110
sims['r801'] = '%s/Eigen/r801_rj95_110_f-'%basedir

#rotated 010
sims['r901'] = '%s/Eigen/r901_rj95_010_f-'%basedir

ray_info={}
ray_info_half = {'ray_type':'obl', 'coord':[0.0,0.505,0.505],'coord2':[0.5,0.505,0.505]}
ray_info_square={'ray_type':'oray','coord':coord,'axis':axis}
ray_info_default=ray_info_square
#this_ray=ray_info_half
#this_ray = ray_info_square
hdx =1./32
two_rays = {'ray_type':'multi', \
        'coord': [ [0.0,0.505,0.505], [hdx,hdx,0.5+hdx]],\
        'coord2':[ [1.0,0.505, 0.505],[np.sqrt(2)+hdx,np.sqrt(2)+hdx, 1.0-hdx]]}
two_other_rays = {'ray_type':'multi', \
        'coord': [ [0.,0.,0.],[0.0,0.505,0.505]],\
        'coord2':[ [1.,1.,1.],[1.0,0.505,0.505]],
        'axis':['x','x']}
ray_x_y = {'ray_type':'multi', \
        'coord': [ [0.0,0.505,0.505], [0.505,0.0,0.505]],\
        'coord2':[ [1.0,0.505, 0.505],[0.505,1.0,0.505]], 'axis':['x','y']}
this_ray = ray_info_square #two_other_rays
sim_list= ['r15','r17']
sim_list=['r15','r15b', 'r17b']
sim_list=['r17b','r17c']
sim_list=['r18b','r19']#,'r18b']
sim_list=['r19b','r19']
sim_list=['r23a','r23b']
sim_list=['r24a','r24b']
#eigen suite
sim_list=['r201','r301']#from enzo, square vs. cos
sim_list=['r201','r401']#square, enzo, rj95
sim_list=['r301','r501']#sin, enzo, rj95
sim_list=['r301','r601']#squares, enzo, rb96
sim_list=['r401','r601']#squares, rj95, rb96
sim_list=['r501','r701']#cos, rj95, rb96
sim_list=['r201','r701']#cos, rj95, rb96
sim_list=['r501','r801']; this_ray=two_rays# cos rj95, 100 vs 110
sim_list=['r501','r901']; this_ray=ray_x_y # cos rj95, 100 vs 110
sim_list=['r202','r402']; this_ray=ray_info_square #left fast
sim_list=['r202','r602']; this_ray=ray_info_square #right fast
sim_list=['r203','r603']; this_ray=ray_info_square #left slow
#sim_list=['r204','r604']; this_ray=ray_info_square #right slow
#sim_list=['r205','r605']; this_ray=ray_info_square #left alf
sim_list=['r206','r606']; this_ray=ray_info_square #left alf
sim_oname = "_%s"*len(sim_list)%tuple(sim_list)
all_ds=[]
frame_list=[50]#,50]
#frame_list=range(0,55,5)
for i, this_frame in enumerate(frame_list):
    frames=[this_frame]

    outname = "Tubec%s_%04d_%s.png"%(sim_oname,this_frame,"_%s"*len(fields)%tuple(fields))

    ds_list = []
    for sim in sim_list:
        ds_list += [yt.load("%s/DD%04d/data%04d"%(sims[sim], frame,frame)) for frame in frames]
    all_ds += ds_list
    print(all_ds)

    g = ds_list[0].index.grids[0]
    data=g
    renorm={'density':1,'pressure':0.6,'TotalEnergy':0.9}
    #renorm = False
    #this_ray = ray_info.get(i,ray_info_default)
    print(this_ray)
    y=tube.tube(ds_list, return_ray=False,  delta=False, fields=fields, renorm=renorm,
            filename = outname,# "Tube_%s%s.png"%(sim,"_%s"*len(fields)%tuple(fields)), 
            #labels=fields,
            plot_args={'marker':'*'},legend=True, labels=sim_list, ylim_method='monotonic', yscale={'density':'linear'}, **this_ray)

#y=tube.tube(all_ds, return_ray=False,  delta=False, fields=fields, renorm=renorm,
#        filename = "TUBE_many",# "Tube_%s%s.png"%(sim,"_%s"*len(fields)%tuple(fields)), 
#        #labels=fields,
#        plot_args={'marker':'*'},legend=True, labels=sim_list+sim_list+sim_list)
