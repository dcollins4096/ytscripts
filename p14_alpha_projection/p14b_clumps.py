

if 'leaf_clumps' not in dir():
    master_clump = Clump(cut_region,"density")
    master_clump.add_validator("min_cells", 20)
    c_min=10
    c_max=cut_region['density'].max()
    step = 2
    find_clumps(master_clump, c_min, c_max, step)
    leaf_clumps = get_lowest_clumps(master_clump) #if both min_cells and grav_bound are used, this is empty.
for ax in [0,1,2]:
    prj=yt.ProjectionPlot(ds,ax,'density',field_parameters=clump_stuff_dict)
    prj.annotate_clumps(leaf_clumps)
    prj.set_cmap('density','gray')
    print prj.save('p14b_t3_clumps_fun')
