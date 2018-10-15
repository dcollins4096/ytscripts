


if 0:
  frame = 50
  basedir = '/scratch1/dcollins/Paper08/B02/512'
  template = basedir+"/RS%04d/restart%04d"
  sim = 'b05'

if 1:
  frame = 0
  basedir = '/scratch1/dcollins/Paper08/B20/512'
  template = basedir+"/RS%04d/restart%04d"
  sim = 'b20'

fname_prefix = "%s_n%04d"%(sim,frame)

if 1:
  ax = 0
  ay = 1
  az = 2
  aa = 0
  d2='xyz'[ay] #cyclic coordinates for Q and U definitions
  d3='xyz'[az]
  ef('p07_functions.py')

if 'ds' not in dir():
  
  """generality comes here.  everything assumes projecting along x"""
  ds = yt.load(template%(frame,frame))
  dproj = ds.proj('density',aa)
  bxproj = ds.proj('Bx',aa)
  byproj = ds.proj('By',aa)
  bzproj = ds.proj('Bz',aa)
  Qproj  = ds.proj("Q",aa)
  Uproj  = ds.proj("U",aa)
  frb_d = dproj.to_frb(1,[512,512])['density']
  frb_bx = bxproj.to_frb(1,[512,512])['Bx']*np.sqrt(4*np.pi)
  frb_by = byproj.to_frb(1,[512,512])['Bx']*np.sqrt(4*np.pi)
  frb_bz = bzproj.to_frb(1,[512,512])['Bx']*np.sqrt(4*np.pi)
  frb_Q = Qproj.to_frb(1,[512,512])['Q']
  frb_U = Uproj.to_frb(1,[512,512])['U']
  bzproj = ds.proj('Bz',aa)

if 1:
  grad_y = np.zeros_like(frb_d)
  grad_z = np.zeros_like(frb_d)
  theta =  np.zeros_like(frb_d)
  grad_y[1:-1,:] = (frb_d[2:,:] - frb_d[:-2,:])/2
  grad_z[:,1:-1] = (frb_d[:,2:] - frb_d[:,:-2])/2

  if 0:
    Phi = np.sqrt(frb_by**2+frb_bz**2)
    BcrossGrad = frb_by*grad_z - frb_bz*grad_y
    BdotGrad =   frb_by*grad_y + frb_bz*grad_z
    theta[1:-1,1:-1] = np.arctan(BcrossGrad[1:-1,1:-1]/BdotGrad[1:-1,1:-1]).v

  if 1:
    Phi_B = 0.5*np.arctan(frb_U/frb_Q)
    Psi = np.zeros_like(frb_d)
    Phi = np.sqrt(frb_U**2+frb_Q**2)
    Psi[1:-1]   = np.arctan(grad_y[1:-1]/grad_z[1:-1])
    theta = Phi_B - Psi.v
    fname_prefix+="_Stokes_"

  Mphi_crit = 0.12/np.sqrt(ds.parameters['GravitationalConstant']/(4*np.pi)) #the 4pis is from Enzo
  MPhi = frb_d/Phi/Mphi_crit
  MPhi_good = (MPhi[1:-1,1:-1].v).flatten()
  theta_good = theta[1:-1,1:-1].flatten()*180/np.pi
  #theta_good[ theta_good<0 ] += 180

if 1:
  """ use good quantities """
  plt.clf()
  nbins = 100
  #hist_val, hist_bins=np.histogram(np.log10(MPhi_good.flatten()), bins=nbins)
  if 0:
      segment_field = np.log10(MPhi_good.flatten()); segment_field_name = "M_phi"
      plt.xlabel(r"$\ln(M/\Phi)/(M/\Phi)_c$")
      plt.ylabel(r"N")
  if 1:
      Sigma = frb_d[1:-1,1:-1].v.flatten()*1.4e22
      segment_field = np.log10(Sigma.flatten()); segment_field_name = "Sigma"
      plt.xlabel(r"$\Sigma [A_{\rm{v}}]$")
      plt.ylabel(r"N")
  hist_val, hist_bins=np.histogram(segment_field, bins=nbins)
  hist_bins_centers=0.5*(hist_bins[:-1]+hist_bins[1:])
  plt.plot(hist_bins_centers, hist_val)

  fname1 = fname_prefix+"hist_"+segment_field_name
  lim_index, cumu_hist = segments(hist_val,5)
  lims = [hist_bins_centers[n] for n in lim_index]
  for n in seg:
      val = hist_bins_centers[n]
      plt.plot([val,val],[0,hist_val.max()])
  plt.savefig(fname1)
  print fname1
  plt.clf()
  hist_args = {'histtype':'step','normed':True,'label':'all','bins':20}
  plt.hist(theta_good,**hist_args)
  hist_args['label'] = "<=%0.2e"%(lims[0])
  plt.hist(theta_good[ segment_field <= lims[0]], **hist_args)
  for n in range(len(lims)-1):
    mask = np.logical_and( segment_field > lims[n], segment_field <= lims[n+1])
    hist_args['label'] = "%0.2d<%s<=%0.2e"%(lims[n],segment_field_name,lims[n+1])
    plt.hist(theta_good[mask], **hist_args)
  hist_args['label'] = ">%0.2e"%(lims[-1])
  plt.hist(theta_good[ segment_field >= lims[-1]], **hist_args)
  plt.legend(loc=3)
  #plt.xlim(-90,120)
  plt.xlabel(r"$\phi$")
  plt.ylabel(r"N")
  fname1 = fname_prefix+"HRO_"+segment_field_name
  plt.savefig(fname1)
  print fname1

# plt.clf()
# Sigma = frb_d[1:-1,1:-1].v.flatten()*1.4e22
# hist_args = {'histtype':'step','normed':True,'label':'all','bins':200}
# plt.clf()
# plt.hist(np.log10(Sigma),**hist_args)
# fname1 = fname_prefix+'Sigma'
# plt.savefig(fname1)
# print fname1


# plt.clf()
# hist_args = {'histtype':'step','normed':True,'label':'all','bins':20}
# plt.hist(theta_good,**hist_args)
# this_field=Sigma
# lims=[10**21.7,10**22.1,10**22.3,10**24.1]
# hist_args['label'] = "<=%0.2e"%lims[0]
# plt.hist(theta_good[ this_field <= lims[0]], **hist_args)
# for n in range(len(lims)-1):
#   mask = np.logical_and( this_field > lims[n], this_field <= lims[n+1])
#   hist_args['label'] = "%0.2e<%s<=%0.2e"%(lims[n],r'$\Sigma$',lims[n+1])
#   plt.hist(theta_good[mask], **hist_args)

# hist_args['label'] = ">%0.2e"%lims[-1]
# plt.hist(theta_good[ this_field <= lims[-1]], **hist_args)

# plt.legend(loc=3)
# #plt.xlim(-90,120)
# plt.xlabel(r"$\phi$")
# plt.ylabel(r"N")
# fname1 = fname_prefix+"phi_sigma"
# plt.savefig(fname1)
# print fname1
# plt.clf()
