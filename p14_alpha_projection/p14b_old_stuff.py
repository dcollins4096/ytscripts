
if 0:
    at.load_clumps('x','cl2d_two',1)
    d=at.leaf_clumps[0]['density']

if 0:
    axis='x'
    clump_prefix='cl2d_two'
    twod_clump = np.transpose(at.get_cl2d(axis,clump_prefix))
    if 'flat_clump' not in dir():
        flat_clump = twod_clump.flatten()
    clump_axis = 0
    ds = yt.load(at.enzo_dataset)
    clump_stuff_tuple = (flat_clump,clump_axis,twod_clump)
    clump_stuff_dict={'clump_mask_stuff':clump_stuff_tuple}
    ds.index.grids[0].set_field_parameter('clump_mask_stuff',clump_stuff_tuple)
    x=ds.index.grids[0]['clump_mask']
    print yt.ProjectionPlot(ds,0,'clump_mask',field_parameters=clump_stuff_dict).save('p14b_clump_mask_c')
    print yt.ProjectionPlot(ds,0,'density').save('p14b_clump_mask_c')

if 0:
    ef('p14_get_cut_region.py')
    ef('p14b_clumps.py')

if 0:
    at.get_cut_region('x','cl2d_two',1)
    #print "should get", 
    print at.cut_region['cell_mass'].sum()
    leafs=at.get_clumps('x','cl2_two',1)

if 1:
    ds=yt.load('/scratch1/dcollins/Paper08//B02/256/RS0030/restart0030')
    ds.index.grids[0].set_field_parameter('bulk_velocity',np.zeros(3))
    sigma2=ds.index.grids[0]['VelocityDispersionSquared']
if 0:
    import mw_stuff
    reload(mw_stuff)
#grav=mw_stuff.run_mw_grav('./mw_test', at.leaf_clumps[0])
    def dump_ascii(mass,x,y,z,filename):
        """yup."""
        fptr = open(filename,'w')
        #print "stuff",filename
        #fptr.write("%15s %15s %15s %15s\n"%("mass","x","y","z"))
        for n in range(len(mass)):
            fptr.write("%0.15e %0.15e %0.15e %0.15e\n"%(mass[n],x[n],y[n],z[n]))
        fptr.close()


if 0:
    import clump_properties
    reload(clump_properties)
    def compute_clump_properties(self):
        pass
        if self.leaf_clumps is None:
            print "run clumps please"
        else:
            self.clump_dir = "%s/clumps_%s"%(self.fits_dir,self.working_clump_prefix)
            clump_stuff_fname = "%s/clump_stuff.pickle"%self.clump_dir
            self.clump_property_list = []

            for n,L in enumerate(self.leaf_clumps):
                self.clump_property_list.append(clump_properties.clump_stuff(L,self.ds, 1),self.clump_index,n)
            fPickle.dump(self.clump_property_list,clump_stuff_fname)
            print "dumped", clump_stuff_fname


compute_clump_properties(at)

