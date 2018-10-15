import tube
reload(tube)
base_dir = '/Users/dcollins/Enzo/enzo-dev/run/MHD/1D/'
bw_ded = '%s/BrioWu-MHD-1D'%base_dir
bw_ct  = '%s/BrioWu-MHD-1D-MHDCT'%base_dir
bw_cbremoval = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/1D/BrioWu-MHD-1D-MHDCT-CenteredBremoval'
bw_ct_fid = '/Users/dcollins/Enzo/dave_vs_enzo_CenteredBremoval/run/MHD/1D/BrioWu-MHD-1D-MHDCT'
frame = 1
fields=['density','pressure','TotalEnergy','x-velocity', 'By' ]#,'x-velocity']
fields += ['ByF','Bz', 'velocity_y','velocity_z']
BWfields = ['density','x-velocity','By']
fields = BWfields
extra=2
frame=11
#set_name = 'ED%02d_%04d/Extra%02d_%04d'%(extra,frame,extra,frame)
#outname = 'ED%02d_%04d.pdf'%(extra,frame)
outname="BW"+"%04d.pdf"%frame
set_name = 'DD%04d/data%04d'%(frame,frame)
#ds_list = [yt.load("%s/%s"%(sim, set_name)) for sim in [bw_ct_fid,bw_cbremoval]]
bw02 = '/scratch1/dcollins/EnzoProjects/E22_DEF_NonCosmo/bw02_ct'
bw01 = '/scratch1/dcollins/EnzoProjects/E22_DEF_NonCosmo/bw01_ded'
outname = 'E22_wtf.pdf'
ds_list = [yt.load("%s/%s"%(sim, set_name)) for sim in [bw01,bw02]]
y=tube.tube(ds_list,  delta=False, fields=fields, renorm=False,filename =outname ,labels=['bw01','bw02'],legend=True)
