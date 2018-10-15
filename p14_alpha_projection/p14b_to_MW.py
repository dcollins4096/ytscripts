

#TO generate text, script for use with Mike Warren's crazy tree code.

def dump_ascii(mass,x,y,z,filename):
    """yup."""
    fptr = open(filename,'w')
    print "stuff",filename
    #fptr.write("%15s %15s %15s %15s\n"%("mass","x","y","z"))
    for n in range(len(mass)):
        fptr.write("%0.15e %0.15e %0.15e %0.15e\n"%(mass[n],x[n],y[n],z[n]))
    fptr.close()

for clump_ind in oober.clumps[frame].keys():
    ef('p14_get_cut_region.py')

    for n_sub,clump in enumerate(oober.clumps[frame][clump_ind].all):
        namedict={}
        #directory = "/lustre/medusa/collins/Paper12/MWgrav/simplegrav"
        directory = "%s/MWgrav"%oober.directory
        name_base = "%s/%s_n%04d_cl%04d_sub%04d"%(directory,oober.name, frame,clump_ind,n_sub)
        #name_base = "%s/TEST"%directory
        ascii_name = name_base+".txt"
        dump_ascii(clump['CellMass'], clump['x'], clump['y'], clump['z'],ascii_name)
        doit_name = name_base + ".doit"
        sdf_name =  name_base + ".sdf"
        namedict['ascii_name'] = ascii_name
        namedict['sdf_name'] = sdf_name
        namedict['ctl_name'] = name_base+".ctl"
        namedict['out_name'] = name_base+".out"
        namedict['name_base'] = name_base
        fptr = open(doit_name,'w')
        fptr.write('grep -v mass %(ascii_name)s | ./txttosdf in=-'%namedict)
        fptr.write(' out=%(sdf_name)s sep=" "  cols=4 desc="float mass, x, y, z"\n'%namedict)
        fptr.write('./nln.x86_64 %(ctl_name)s > %(out_name)s\n'%namedict)
        fptr.close()
        fptr = open(namedict['ctl_name'],'w')
        fptr.write('# SDF\n')
        fptr.write('char datafile[] = "%(sdf_name)s";\n'%namedict)
        fptr.write('char outfile[] = "%(name_base)s";\n'%namedict)
        fptr.write('int save_first = 1;\n')
        fptr.write('float dt = 0.01;\n')
        fptr.write('int do_Arel = 1;\n')
        fptr.write('int do_DL = 1;\n')
        fptr.write('float epsilon = 1e-4;\n')
        fptr.write('float errtol = 0.0001;\n')
        fptr.write('float frac_tol = 0.0001;\n')
        fptr.write('int nsteps = 0;\n')
        fptr.write('int Msg_memfile = 65536;\n')
        fptr.write('# SDF-EOH\n')
        fptr.close()




        
