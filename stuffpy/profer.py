#!/usr/bin/env python
fname_list = []
plt.clf()
clist = ['r','g','b']

class the_stuff():
    def __init__(self,name,fname,color):
       self.simname = fname.split('/')[-1]
       self.color=color
       self.fname=fname
       self.name=name
       self.fraction=[]
       self.cuml=[]
       self.self=[] #the one I want
       self.calls=[]
       self.self_ms_per_call=[]
       self.total_ms_per_call=[]
       self.name=[]
       self.verbose=False
    def scrub(self):
        fptr = open(self.fname)
        lines = fptr.readlines()
        for line in lines[5:]:
            if self.verbose:
                print "==================="
            self.fraction.append(float(line[0:7]))
            if self.verbose:
                print self.fraction[-1]
            self.cuml.append(float(line[7:17]))
            if self.verbose:
                print self.cuml[-1]
            self.self.append(float(line[17:26]))
            if self.verbose:
                print self.self[-1]
            try:
                self.calls.append(int(line[26:35]))
            except:
                self.calls.append(-1)
            if self.verbose:
                print self.calls[-1]

            try:
                self.self_ms_per_call.append(float(line[35:44]))
            except:
                self.self_ms_per_call.append(-1)
            if self.verbose:
                print self.self_ms_per_call[-1]
            try:
                self.total_ms_per_call.append(float(line[44:53]))
            except:
                self.total_ms_per_call.append(-1)
            if self.verbose:
                print self.total_ms_per_call[-1]
            this_name_full = line[54:]
            if self.verbose:
                print this_name_full
            this_name = this_name_full.split('(')[0]
            if this_name[-1] == "\n":
                this_name = this_name[:-1]
            if self.verbose:
                print this_name
            self.name.append(this_name)
            if self.self[-1] < 0.009:
                return
    def plot(self,pd, names=True,cuml=False):
        if cuml:
            to_plot=np.cumsum(self.self)
        else:
            to_plot = self.self
        pd.plot(to_plot,marker="+", label=self.simname,color=self.color)
        if names:
            for n in range(len(self.self)):
                plt.text(n,to_plot[n],self.name[n],fontsize=5,color=self.color)

if 0:
    gpd = the_stuff('gpd','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b11_Timing_Ded/gprofD2','r')
    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b10_Etiming_F/gprofF2','g')
    gpc = the_stuff('gpc','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b09_Etiming_C/gprofC2','b')
if 1:
#    gpc = the_stuff('gpd','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b12_old_C/gprof_old_c','r')
    gpc = the_stuff('gpd','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b12_old_C/gprof_old_c_long','r')
#    gpc = the_stuff('gpd','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b12_old_C/gprof_old_c_2','r')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b13_old_F/gprof_old_f','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_prim1','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_press','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_limiter1','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_limiter2','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_limiter3','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_new_timestep','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_nt_2','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_center','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_cons1','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_center_fortran','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_limiter_method0','g')
    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_long','g')
#    gpf = the_stuff('gpf','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b15_old_F_opt/gprof_old_f_limiter4','g') #hack it out.
#    gpd = the_stuff('gpc','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b14_old_D/gprof_old_d','b')
#    gpd = the_stuff('gpc','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b14_old_D/gprof_old_d_2','b')
    gpd = the_stuff('gpc','/scratch1/dcollins/EnzoProjects/E27_MHDPerformance/b14_old_D/gprof_old_d_long','b')
gpd.scrub()
gpc.scrub()
gpf.scrub()
plt.clf()
names=True
nplot=15

if nplot>0:
    if 0:
        cuml=True
        outname = 'cuml%02d.pdf'%nplot
        loc=4
    else:
        cuml=False
        outname = 'straigt%02d.pdf'%nplot
        loc=1
    gpd.plot(plt,names,cuml)
    gpc.plot(plt,names,cuml)
    gpf.plot(plt,names,cuml)
    plt.legend(loc=loc)
    plt.savefig(outname)
    print outname
ef('proferb.py')
print "gpd", sum(gpd.self)
print "gpc", sum(gpc.self)
print "gpf", sum(gpf.self)
