
from volavg import *


class ZeroField():
    def __init__(self,  refine_by=1):
        self.fields_to_replace = ['Bx','By','Bz','BxF','ByF','BzF','TotalEnergy']
    def NewFunction(self,target_grid, field):
        return target_grid[field][:]*0

class FillerFunction():
  #kludge refine by, domain dimension
  def __init__(self,pf,refine_by=1):
      self.cg = pf.h.covering_grid(0,[0.0]*3,pf.domain_dimensions)
      self.small={}
      self.extra_dims = {'BxF':nar([1,0,0]),'ByF':nar([0,1,0]),'BzF':nar([0,0,1])}
      self.fields_to_replace=[]
      fine_density = self.cg['Density']
      coarse_density = volavg(fine_density, rank=3, refine_by = refine_by)
      for in_field_name in ['BxF', 'ByF', 'BzF','density',
                           'x-velocity','y-velocity','z-velocity',
                           'Bx','By','Bz']: #['Density','x-velocity','y-velocity','z-velocity','BxF','ByF','BzF']:
          infield = self.cg[in_field_name]
          dims = infield.shape
          #out_field_name = {'MagneticField_F_1':'BxF','MagneticField_F_2':'ByF','MagneticField_F_3':'BzF',
          #                'Density':'Density','x-velocity':'x-velocity','y-velocity':'y-velocity','z-velocity':'z-velocity'}[in_field_name]
          out_field_name = in_field_name
          self.fields_to_replace.append(out_field_name)
          if out_field_name is 'BxF':
            field_to_store = na.ones( [(dims[0])/refine_by+1, dims[1]/refine_by, dims[2]/refine_by])
            for i in range(dims[0]/refine_by):
                field_to_store[i,:,:] = volavg(infield[refine_by*i,:,:], rank=2, refine_by=refine_by)
            field_to_store[dims[0]/refine_by,:,:] = field_to_store[0,:,:]
          elif out_field_name is 'ByF':
            field_to_store = na.ones( [(dims[0])/refine_by, (dims[1])/refine_by+1, dims[2]/refine_by])
            for j in range(dims[1]/refine_by):
                field_to_store[:,j,:] = volavg(infield[:,refine_by*j,:], rank=2, refine_by=refine_by)
            field_to_store[:,dims[1]/refine_by,:] = field_to_store[:,0,:]
          elif out_field_name is 'BzF':
            field_to_store = na.ones( [(dims[0])/refine_by, dims[1]/refine_by, (dims[2])/refine_by+1])
            for k in range(dims[2]/refine_by):
                field_to_store[:,:,k] = volavg(infield[:,:,refine_by*k], rank=2, refine_by=refine_by)
            field_to_store[:,:,dims[2]/refine_by] = field_to_store[:,:,0]
          elif out_field_name is 'Density':
            field_to_store = coarse_density
          else:
              field_to_store = volavg(fine_density*infield,rank=3,refine_by=refine_by)/coarse_density
          self.small[out_field_name] = field_to_store.swapaxes(0,2)
  def NewFunction(self,target_grid, field):
      s= target_grid.get_global_startindex()
      e= s + target_grid.ActiveDimensions + self.extra_dims.get(field,na.zeros(3))
      #output = self.small[field][ s[0]:e[0], s[1]:e[1], s[2]:e[2]]
      output = self.small[field][ s[2]:e[2], s[1]:e[1], s[0]:e[0]]  #because swapaxes.
      return output

def refile(pf,Filler):
    """For each grid in *pf*, use *Filler* to replace the data with whatever *Filler* returns"""
    print Filler.fields_to_replace
    print "NEW THING"
    for g in pf.index.grids:
        print g
        fptr = h5py.File(g.filename,'r+')
        print g.filename
        grd = fptr['Grid%08d'%g.id]
        for field in grd.keys():
            if field in Filler.fields_to_replace:
                #output = Filler.NewFunction(g,field)
                if field == 'TotalEnergy':
                    output = grd[field][:] -  0.5*(1e-15)**2/grd['Density'][:]
                else:
                    output = grd[field][:]*0
                if grd[field].shape != output.shape:
                    print "CLOWN falure"
                else:#
                    print "i am so smart", field
                    del grd[field]
                    grd.create_dataset(field,data=output)
                
        fptr.close()

g19e_pf = yt.load('/mnt/c/scratch/sciteam/dcollins/Paper33_Galaxy/g19e_R128L4_h6_b0/DD0000/DD0000')
filler=ZeroField() 
refile(g19e_pf,filler)
if 'SmallFive' not in dir() and False:
    #Filler Function takes pf_large and prepares the small grids.
    #refile takes the FillerFunction and refills pf_target
    pf_large = yt.load('/scratch1/dcollins/Paper08/B2/512/RS0000/restart0000')
    SmallFive = FillerFunction(pf_large,refine_by=4)
    pf_target = yt.load('/scratch1/dcollins/Paper08/B2/128/DD0000/data0000')
#refile(pf_target,SmallFive)
#pf_ok = load('/scratch1/dcollins/Paper19/B02/128c/DD0001/data0001')
#ProjectionPlot(pf_ok,0,'MagneticEnergy').save('d1b')



