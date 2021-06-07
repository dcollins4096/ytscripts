from go import *
from collections import defaultdict
import sim_colors
reload(sim_colors)
output_dir = os.environ['HOME']+'/PigPen/'

sims=['1_1','1_half','2_2','3_1','3_half','half_2', 
      '1_2','2_1','2_half','3_2','half_1','half_half']
sims=sims

field_list=[ 'Bx_fluct','By_fluct','Bz_fluct']

PDFs = {}
Bins = {}
for sim in sims:
    PDFs[sim]=defaultdict(list)
    Bins[sim]=defaultdict(list)

    for field in field_list:
        file_list = glob.glob("ProfileFiles/p49_%s_????_%s*h5"%(sim,field))

        for fname in file_list:
            bins,pdf = dpy(fname, [field,'cell_volume'])
            PDFs[sim][field].append(pdf)
            Bins[sim][field].append(bins)
        PDFs[sim][field] = nar(PDFs[sim][field])
        PDFs[sim][field] = PDFs[sim][field].mean(axis=0)
        Bins[sim][field] = nar(Bins[sim][field])
        Bins[sim][field] = Bins[sim][field].mean(axis=0)

fig, axgrid = plt.subplots(2,2)
axlist=axgrid.flatten()


for sim in sims:
    for nf, field in enumerate(field_list):
        ls = sim_colors.linestyle[sim]
        c = sim_colors.color[sim]
        bins = 0.5*(Bins[sim][field][1:] + Bins[sim][field][:-1])
        db   =     (Bins[sim][field][1:] - Bins[sim][field][:-1])
        axlist[nf].plot( bins, PDFs[sim][field]/db, linestyle=ls, c=c)
        axbonk(axlist[nf],xlabel=field,ylabel='V',xscale='linear',yscale='log')
        #axlist[nf].set_xscale('symlog',linthresh=1e-5)
        
fig.savefig('%s/all_pdf.pdf'%output_dir)



