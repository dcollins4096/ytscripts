"""
Code to facilitate use of Mike Warren's sweet and very fast potential solver.
"""
import os
import subprocess
def mw_grav_write(name_base, clump, execpath="./"):
    namedict={}
    ascii_name = name_base+".txt"
    dump_ascii(clump['cell_mass'].in_units('code_mass'), clump['x'].in_units('code_length'), clump['y'].in_units('code_length'), clump['z'].in_units('code_length'),ascii_name)
    doit_name = name_base + ".doit"
    sdf_name =  name_base + ".sdf"
    namedict['ascii_name'] = ascii_name
    namedict['sdf_name'] = sdf_name
    namedict['ctl_name'] = name_base+".ctl"
    namedict['out_name'] = name_base+".out"
    namedict['name_base'] = name_base
    namedict['execpath'] = execpath
    fptr = open(doit_name,'w')
    fptr.write('#!/bin/tcsh\n')
    fptr.write('grep -v mass %(ascii_name)s | %(execpath)s/txttosdf in=-'%namedict)
    fptr.write(' out=%(sdf_name)s sep=" "  cols=4 desc="float mass, x, y, z"\n'%namedict)
    fptr.write('%(execpath)s/nln.x86_64 %(ctl_name)s > %(out_name)s\n'%namedict)
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

def run_mw_grav(name_base,clump):
    doit_name = name_base + ".doit"
    execpath = "/home/dcollins4096/Scripts"
    mw_grav_write(name_base,clump, execpath=execpath)
    os.chmod(doit_name,0744)
    subprocess.call(doit_name)
    outfile = name_base+".out"
    fptr = open(outfile,'r')
    for line in fptr:
        if line.startswith('ke'):
            mw_grav = float(line.split(" ")[5][:-1])
            break
    fptr.close()
    return mw_grav
    
    


def dump_ascii(mass,x,y,z,filename):
    """yup."""
    fptr = open(filename,'w')
    #print "stuff",filename
    #fptr.write("%15s %15s %15s %15s\n"%("mass","x","y","z"))
    for n in range(len(mass)):
        fptr.write("%0.15e %0.15e %0.15e %0.15e\n"%(mass[n],x[n],y[n],z[n]))
    fptr.close()


