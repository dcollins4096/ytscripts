reload(p49_plot_tools)

fb=p49_eigen.get_ffts(tsb.cubes)
ff=p49_eigen.get_ffts(tsf.cubes)
fs=p49_eigen.get_ffts(tss.cubes)


if 0:
    field='d'
    print("== fb ==")
    print_nz_2(fb[field])
    print("== ff ==")
    print_nz_2(ff[field])
    print("== fs ==")
    print_nz_2(fs[field])
if 1:
    sb = p49_plot_tools.chomp(out_dir_b)
    sf = p49_plot_tools.chomp(out_dir_f)
    ss = p49_plot_tools.chomp(out_dir_s)
    out_b="./AAA/test_b_"
    out_s="./AAA/test_s_"
    out_f="./AAA/test_f_"
    analysis={'print_wave':True,
              'plot_fields':0,
              'k_mag':0,
              'k_proj':0}
    print("===== both =====")
    p49_plot_tools.do_stuff(stuff=sb,outdir=out_b,**analysis)
    print("===== fast =====")
    p49_plot_tools.do_stuff(stuff=sf,outdir=out_s,**analysis)
    print("===== slow =====")
    p49_plot_tools.do_stuff(stuff=ss,outdir=out_f,**analysis)

if 0:
    field='d'
    print("== fb f- ==")
    print_nz_2(tsb.hat_box['f-'][field])
    print("== fb s+ ==")
    print_nz_2(tsb.hat_box['s+'][field])
    print("== ff f- ==")
    print_nz_2(tsf.hat_box['f-'][field])
    print("== fs s+ ==")
    print_nz_2(tss.hat_box['s+'][field])
#   print("== ff ==")
#   print_nz_2(ff[field])
#   print("== fs ==")
#   print_nz_2(fs[field])

if 0:
    sb={'cubes':tsb.cubes,'means':tsb.quan}
    sf={'cubes':tsf.cubes,'means':tsf.quan}
    ss={'cubes':tss.cubes,'means':tss.quan}
    p49_plot_tools.print_wave_content(stuff=sb)
    p49_plot_tools.print_wave_content(stuff=sf)
    p49_plot_tools.print_wave_content(stuff=ss)

if 0:
    sb = p49_plot_tools.chomp(out_dir_b)

if 0:
    sb['means']=tsb.quan

    analysis={'print_wave':True,
              'plot_fields':0,
              'k_mag':0,
              'k_func':maxis,
              'k_proj':0}
    p49_plot_tools.do_stuff(stuff=sb,outdir="./AAA/redo_b_",**analysis)
field_list = ['d','px','py','pz','hx','hy','hz','p','vx','vy','vz']
field_list = ['d','hx','hy','hz','p','vx','vy','vz']
#for field in field_list:
#    print("%5s %0.2e"%(field,np.var(sn['cubes'][field]-tsb.cubes[field])))
for x in []:# 'xyz':
    vn = 'v'+x
    v1=sb['cubes']['v'+x]
    p1=sb['cubes']['p'+x]
    d1=sb['cubes']['d']
    v2 = tsb.cubes['v'+x]
    p2 = tsb.cubes['p'+x]
    d2 = tsb.cubes['d']
    vf = tsf.cubes['v'+x]
    pf = tsf.cubes['p'+x]
    df = tsf.cubes['d']
    vs = tss.cubes['v'+x]
    ps = tss.cubes['p'+x]
    ds = tss.cubes['d']
    print("%5s %0.2e"%(x,np.var((1.-d2)-(1.-df)-(1.-ds))))
