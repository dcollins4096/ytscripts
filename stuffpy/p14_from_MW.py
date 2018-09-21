
for clump_ind in oober.clumps[frame].keys():
    for n_sub in range(1,len(oober.clumps[frame][clump_ind].all)):
        clump = oober.clumps[frame][clump_ind].all[n_sub]
        #directory = '/lustre/medusa/collins/Paper12/MWgrav/simplegrav'
        directory = "%s/MWgrav"%oober.directory
        name_base = "%s/%s_n%04d_cl%04d_sub%04d"%(directory,oober.name, frame,clump_ind,n_sub)
        print "%s_n%04d_cl%04d_sub%04d"%(oober.name, frame,clump_ind,n_sub)
        if not hasattr(clump,'stuff') or clump.stuff is None:
            print "gotta make some stuff."
            clump.stuff= clump_stuff.clump_stuff_code(clump,10)

        outname = name_base+".out"
        G = clump.data.pf['GravitationalConstant']/(4*na.pi)
        fptr = open(outname,'r')
        for line in fptr:
            if line.startswith('ke'):
                mw_grav = -G*float(line.split(" ")[5][:-1])
                break
        fptr.close()
        clump.stuff.Energies['MWgrav'] = mw_grav
        if clump.stuff.Energies['Gravitational'] > 0:
            print relerr(clump.stuff.Energies['MWgrav'],clump.stuff.Energies['Gravitational']) 
        else:
            print "no old", clump.stuff.Energies['MWgrav']

        clump.stuff.Ratios_mw = {}
        ge = mw_grav
        if ge < 0:
            print "WTF?"
            pdb.set_trace()
        ke = clump.stuff.Energies['Kinetic']
        be = clump.stuff.Energies['Magnetic']
        thermal = clump.stuff.Energies['Thermal']
        if ke+be+thermal > 0.0:
            clump.stuff.Ratios_mw['ge_mw/all']=ge/(ke+be+thermal)
        if ke+be > 0.0:
            clump.stuff.Ratios_mw['ge_mw/ke+be']=ge/(ke+be)
        if ke > 0.0:
            clump.stuff.Ratios_mw['ge_mw/ke']=ge/ke
        if be > 0.0:
            clump.stuff.Ratios_mw['ge_mw/be']=ge/be
            clump.stuff.Ratios_mw['ke/be']=ke/be
            clump.stuff.Ratios_mw['te/be'] = thermal/be
        if thermal > 0.0:
            clump.stuff.Ratios_mw[ 'ge_mw/thermal']=ge/thermal
        if ke+ge > 0.0:
            clump.stuff.Ratios_mw['ge_mw/(ke+ge)'] = ge/(ke+ge)
        if clump.stuff.Energies['Rotational'] > 0.0:
            clump.stuff.Ratios_mw['ge_mw/re'] = clump.stuff.Energies['MWgrav']/clump.stuff.Energies['Rotational']

    if 0:
        new_clump_name = "%s/%s_%04d/%s_%04d"%(oober.directory,'clumps_MW',frame, 'clump', clump_ind)
        oober.clumps[frame][clump_ind].save(new_clump_name)
    else:
        print "NOT SAVING"



if 0:
    """stick this in the fptr bit to get timing."""
    if line.startswith('Step'):
        counter = 0
        for i in line.split(" "):
            if i != "":
                counter += 1
            if counter == 3:
                break
        optr.write( " %s\n"%i)
