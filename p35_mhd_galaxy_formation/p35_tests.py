

base_dir = "/mnt/c/scratch/sciteam/dcollins/Paper35_cosmology"

all_sims = glob.glob("%s/*"%(base_dir))
hard_sim = ['/mnt/c/scratch/sciteam/dcollins/Paper35_cosmology/b04_small']
hard_sim = [base_dir+'/c04_small_medium_mhd']
for sim in hard_sim:
    red10 = sim+"/RD%04d/RD%04d"%(10,10)
    if len(glob.glob(red10)) > 0:
        ds = yt.load(red10)
        proj = yt.ProjectionPlot(ds,0,'density')
        proj.annotate_grids()
        print proj.save(sim)

