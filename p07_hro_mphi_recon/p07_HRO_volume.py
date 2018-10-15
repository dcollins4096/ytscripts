
if 1:
  frame = 10
  basedir = '/scratch1/dcollins/Paper08/B02/512'
  template = basedir+"/RS%04d/restart%04d"
  sim = 'b02'

if 0:
  frame = 50
  basedir = '/scratch1/dcollins/Paper08/B20/512'
  template = basedir+"/RS%04d/restart%04d"
  sim = 'b20'

fname_prefix = "%s_n%04d"%(sim,frame)
if 'ds' not in dir():
  ds = yt.load(template%(frame,frame))
if 'cg' not in dir():
  cg = ds.covering_grid(0,[0]*3,[512]*3)
if 1:
  den = cg['density'].v*1000
  grad_x = np.zeros_like(den)
  grad_y = np.zeros_like(den)
  grad_z = np.zeros_like(den)
  grad_x[1:-1,:,:] = den[2:,:,:]-den[:-2,:,:]
  grad_y[:,1:-1,:] = den[:,2:,:]-den[:,:-2,:]
  grad_z[:,:,1:-1] = den[:,:,2:]-den[:,:,:-2]
  bx = cg['Bx'].in_units('code_magnetic').v*np.sqrt(4*np.pi)
  by = cg['By'].in_units('code_magnetic').v*np.sqrt(4*np.pi)
  bz = cg['Bz'].in_units('code_magnetic').v*np.sqrt(4*np.pi)
  btot = np.sqrt(bx*bx+by*by+bz*bz)
  Mphi_crit = 0.12/np.sqrt(ds.parameters['GravitationalConstant']/(4*np.pi)) #the 4pis is from Enzo
  MPhi = (den*1./512)/btot/Mphi_crit

  BcrossGrad = (by*grad_z-bz*grad_y)**2+\
               (bz*grad_x-bx*grad_z)**2+\
               (bx*grad_y-by*grad_x)**2
  BcrossGrad = np.sqrt(BcrossGrad)
  BdotGrad = bx*grad_x + by*grad_y+bz*grad_z
  theta = np.zeros_like(by)
  theta[1:-1,1:-1,1:-1] = np.arctan(BcrossGrad[1:-1,1:-1,1:-1]/BdotGrad[1:-1,1:-1,1:-1])
  MPhi_good = MPhi[1:-1,1:-1,1:-1].flatten()
  theta_good = theta[1:-1,1:-1,1:-1].flatten()
# theta_good *= 180./np.pi
# theta_good[theta_good<0] += 180

if 1:
  plt.clf()
  plt.hist(np.log10(MPhi_good.flatten()),bins=100)
  plt.xlabel(r"$\ln(M/\Phi)/(M/\Phi)_c$")
  plt.ylabel(r"N_{vol}")
  fname1 = fname_prefix+"Mphi_vol"
  plt.savefig(fname1)
  print fname1
            
  plt.clf()
  plt.hist(np.log10(den.flatten()),bins=100)
  plt.xlabel(r"$\ln(n)$")
  plt.ylabel(r"N_{vol}")
  fname1 = fname_prefix+"n_vol"
  plt.savefig(fname1)
  print fname1
  plt.clf()

  plt.clf()
  plt.hist(theta_good,bins=100)
  plt.xlabel(r"$\phi$")
  plt.ylabel(r"N_{vol}")
  fname1 = fname_prefix+"phi_vol"
  plt.savefig(fname1)
  print fname1

if 1:
  plt.clf()
  this_field=den[1:-1,1:-1,1:-1].flatten()
  this_x = np.cos(theta_good)
  #this_x = theta_good
  hist_args = {'histtype':'step','normed':True,'label':'all','bins':20}
  lims=[10**2,10**3,10**4,10**5]
  hist_args['label'] = "<=%0.0e"%lims[0]
  plt.hist(this_x[ this_field <= lims[0]], **hist_args)
  print (this_field <= lims[0]).shape
  for n in range(len(lims)-1):
    mask = np.logical_and( this_field > lims[n], this_field <= lims[n+1])
    hist_args['label'] = "%0.0e<%s<=%0.0e"%(lims[n],r'$n$',lims[n+1])
    plt.hist(this_x[mask], **hist_args)
    print this_x[mask].shape, this_x[mask].min(), this_x[mask].max()
  hist_args['label'] = ">%0.0e"%lims[-1]
  plt.hist(this_x[ this_field >= lims[-1]], **hist_args)
  plt.legend(loc=1)
  #plt.xlim(-90,120)
  plt.xlabel(r"$\phi$")
  plt.ylabel(r"N")
  fname1 = fname_prefix+"phi_cut"
  plt.savefig(fname1)
  print fname1
  plt.clf()

if 1:
  plt.clf()
  this_field=MPhi_good
  this_x = np.cos(theta_good)
  #this_x = theta_good
  hist_args = {'histtype':'step','normed':True,'label':'all','bins':20}
  lims=[1,5,10,100]
  hist_args['label'] = "<=%0.0e"%lims[0]
  plt.hist(this_x[ this_field <= lims[0]], **hist_args)
  print (this_field <= lims[0]).shape
  for n in range(len(lims)-1):
    mask = np.logical_and( this_field > lims[n], this_field <= lims[n+1])
    hist_args['label'] = "%0.0e<%s<=%0.0e"%(lims[n],r'$M/\Phi$',lims[n+1])
    plt.hist(this_x[mask], **hist_args)
    print this_x[mask].shape, this_x[mask].min(), this_x[mask].max()
  hist_args['label'] = ">%0.0e"%lims[-1]
  plt.hist(this_x[ this_field >= lims[-1]], **hist_args)
  plt.legend(loc=1)
  #plt.xlim(-90,120)
  plt.xlabel(r"$\phi$")
  plt.ylabel(r"N")
  fname1 = fname_prefix+"phi_cut_MPhi"
  plt.savefig(fname1)
  print fname1
  plt.clf()
