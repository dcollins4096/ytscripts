
import re
def vis_file(filename):
  fptr = open(filename,'r')
  print "started", filename
  #levelstep = re.compile(r'Level[0]: dt = 1.39534e-05\(This So Far: 1.39534e-05 Above: 1.39534e-05 Frac 1\)')
  string = r'^Level\[(\d)\]: dt = (\d\.\d+....)\(This So Far: (\d\.\d+....) Above: (\d\.\d+....).*'
  #string = r'^Level\[(\d)\]: dt = (\d\.\d+.......)\s*(\d.\d\d\d\d\d\d\d)\s*\((\d.\d\d\d\d\d\d\d)/(\d.\d\d\d\d\d\d\d)\)'
  #string = r'^Level\[(\d)\]: dt = (\d\.\d+)\s*(\d.\d+)\s*\((\d\.\d+)/(\d.\d\d\d\d\d\d\d)\)'
  #string = r'^Level\[(\d)\]: dt = (\d\.\d+)\s*(\d.\d+)\s*\(([^/]*)/([^/]*)\)'
  string = r'^Level\[(\d)\]: dt = ([^\s]+)\s*([^\s]+)\s*\(([^\s]+)/([^\s]+)\)'
#                     1              2                   3                  !   4                    5
  #string = r'^Level\[(\d)\]: dt = (\d\.\d+.......)   \d.\d\d\d\d\d\d\d \(\d.\d\d\d\d\d\d\d/\d.\d\d\d\d\d\d\d\)'
  """
  Level[0]: dt = 0.0183712
  Level[1]: dt = 0.0183712  0.0183712 (0.0183712/0.0183712)
  """
  levelstep = re.compile(string)
  all_the_things=[]
  level_time = {}
  all_time=[]
  all_level=[]
  level_time_rel = {}
  level_ind = {0:[],1:[],2:[],3:[],4:[],5:[]}
  color_list=[]
  counter = 0
  for line in fptr:
    #print line[:-1]
    if counter < -1:
      break
    match = levelstep.match(line)
    if match is None:
        if line.startswith('Level'):
            print line
    else:
      #print ".",
      if counter > 10000:
        break
    
      level = int(match.group(1))
      dt = float(match.group(2))
      SoFar = -1 #float(match.group(3))
      Above = float(match.group(5))
      #print level,dt, SoFar, Above, dt/Above
      if level_time.has_key(level):
          level_time[level] +=  dt
      elif level_time.has_key(level-1):
          parent_start_time = level_time[level-1] #- Above
          level_time[level] = parent_start_time
          #print "NL1:", parent_start_time, Above, level_time[level-1]
      else:
          level_time[level] = 0
          print "NL0: 0"
      list_to_plot = [counter, level, level_time[level], dt,dt/Above]
      k=[0.0,0.0,0.0,1.0]
      r=[1.0,0.0,0.0,1.0]
      g=[0.0,1.0,0.0,1.0]
      b=[0.0,0.0,1.0,1.0]
      c=[0.5,0.5,0.0,1.0]
      m=[1.0,0.0,1.0,1.0]
      y=[1.0,1.0,0.0,1.0]
      gray=[0.5,0.5,0.5,1]
      color_list.append([k,r,g,b,c,m,y,gray,gray][level])
      #print list_to_plot
      all_the_things.append(list_to_plot)
      #level_ind[level].getappend(counter)
      counter += 1
      #print "...",counter
      #print all_the_things
      #print "."
      #print list_to_plot, len(all_the_things), len(all_the_things[-1])
      x = nar(all_the_things)


  all_the_things = nar(all_the_things)
  for l in level_ind.keys():
    level_ind[l] = nar(level_ind[l])
  plt.clf()
  #plt.plot( all_the_things[:,2])
  return all_the_things,color_list, level_ind

if 1:
    
    #fname_ak04 = '/Users/dcollins/RESEARCH2/Paper35_CosmologyMHD/2016-08-18-dt/5311520.bw.OU'; oname = 'ak04'
    #fname_ak04 = '/Users/dcollins/RESEARCH2/Paper35_CosmologyMHD/2016-08-18-dt/ak01.o5287803'; oname = 'ak01' #woah, so long.
    fname_ak04 = '/Users/dcollins/RESEARCH2/Paper35_CosmologyMHD/2016-08-18-dt/ak01.o.Level'; oname = 'ak01'
    all_ak04, color_ak04, ind_ak04 = vis_file(fname_ak04)
    color_ak04=nar(color_ak04)
    print all_ak04
if 1:
    plt.clf()
    mslice = slice(None)
    
    if 1:
        TheX = all_ak04[mslice,0]; TheXLabel = 'Nevent'
        TheY = all_ak04[mslice,3]; TheYLabel = 'dt'
        output = '%s_%s_%s.pdf'%(oname,'Nevent','dt')
    plt.plot(TheX,TheY ,c=[0.5]*3)
    plt.scatter(TheX, TheY,marker='*',c=color_ak04,linewidths=0, s=100)
    plt.xlabel(TheXLabel); plt.ylabel(TheYLabel);
    plt.yscale('log')
    plt.savefig(output)
    print output
if 0:
    """old stuff with P27"""
#filename = '/Users/dccollins/ThisScratch/Paper12/B20/256_j4/Meta/B20_j4.o878433'


    if 0:
      all,color,all_level_ind=vis_file(filename)
      cycles=all[:,0]
      dt = all[:,3]
      ratio = all[:,4]

    if 0:
      """The first root grid timestep"""
      TheSlice = slice(0,105)

    CompleteSteps = 1078
    if 0:
      """11 complete timesteps."""
      TheSlice = slice(0,CompleteSteps)

    if 0:
      """level_ind consistent with slices"""
      level_ind={}
      for L in range(5):
        level_ind[L] = all_level_ind[L][ all_level_ind[L] < CompleteSteps]


#TheY = all[TheSlice,1]
#TheY = all[TheSlice,4]
    TheX = None
    TheX_Lim = None
    TheY_Lim = None
    TheX_Scale = 'linear'
    TheY_Scale = 'linear'

    if 0:
      
      plt.clf()
      dt_0 = dt[level_ind[0]]
      plt.plot(dt_0)
      for l in range(5):
        print "total time, level %d",l, dt[level_ind[l]].sum()

    if 0:
      for l in range(5):
        print "nsteps, level %d", l, len(level_ind[l])
      print "nsteps, total", sum( [len(level_ind[l]) for l in range(5)] )

    if 0:
      """assuming that the first step is representative, which is is not,
      which level dominates the timestep?
      If we enforce dt_l = """
      total_time = dt[level_ind[0]].sum()
      first_steps = nar( [ dt[level_ind[L][0]] for L in range(5) ])
      print "Step length", "%0.2e, "*len(first_steps)%tuple(first_steps)
      print "Step length, rel  0", "%0.2e, "*len(first_steps)%tuple(first_steps/first_steps[0])
      print "Step length, rel -1 %8s,"%('--'), "%0.2e, "*(len(first_steps)-1)%tuple(first_steps[1:]/first_steps[0:-1])
      ideal_root = first_steps * 2**(na.arange(5))
      print "Ideal Root", "%0.2e, "*len(ideal_root)%tuple(ideal_root)
      print "Ideal Root", "%0.2e, "*len(ideal_root)%tuple(ideal_root/ideal_root[0])
      min_ideal = min(ideal_root)
      print "Min ideal %0.2e / %0.2e = %0.2e"%(min_ideal,first_steps[0],min_ideal/first_steps[0])
      print "Total %0.2e"%total_time
      n_ideal_steps = 0
      n_steps_total = sum( [len(level_ind[l]) for l in range(5)] )
      for L in range(5):
        print "nsteps at min_ideal %d %e"%(L,total_time/(min_ideal*0.5**(L)))
        n_ideal_steps += total_time/(min_ideal*0.5**(L))
      print "nsteps, ideal %e, reduced by %e"%(n_ideal_steps, n_ideal_steps/n_steps_total)
      

    if 0:
      """legend"""
      k=[1.0,1.0,1.0,1.0]
      r=[1.0,0.0,0.0,1.0]
      g=[0.0,1.0,0.0,1.0]
      b=[0.0,0.0,1.0,1.0]
      c=[0.5,0.5,0.0,1.0]
      circle_list = [plt.scatter([0],[0],c=d) for d in [k,r,g,b,c]]
      plt.clf()

    if 0:
      """x info"""
      TheX_Lim = -1,128
      TheX = cycles[TheSlice]
      TheX_Label = 'Call'
      
    if 1:
      """dt, linear"""
      TheY = dt[TheSlice]
      TheY_Scale = 'linear'
      TheY_Lim = 1e-7,1.5e-5
      TheY_Label = 'dt'

    if 0:
      """dt, log"""
      TheY = dt[TheSlice]
      TheY_Scale = 'log'
      TheY_Lim = 1e-7,1.5e-5
      TheY_Label = 'dt'

    if 0:
      TheC = nar(color[TheSlice])
      OughtToBe=TheY[0]*0.5**(na.arange(5))
      plt.scatter([TheX[0]]*len(OughtToBe),OughtToBe, color=nar([k,r,g,b,c]))

    if 0:
      TheY = ratio[TheSlice]
      TheY_Scale = 'linear'
      TheY_Lim = None
      TheY_Label='ratio'

    if 0:
        if TheX is not None:
          plt.yscale(TheY_Scale)
          #plt.ylim(1e-7,1e-4)
          plt.ylim(TheY_Lim)
          plt.xlim(TheX_Lim)
          plt.scatter(TheX, TheY, c=TheC)
          plt.legend(circle_list,range(5))
          plt.ylabel(TheY_Label)
          plt.xlabel(TheX_Label)
