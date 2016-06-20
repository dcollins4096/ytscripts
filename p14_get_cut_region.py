if 'clump_mask_defined' not in dir():
    ef('p14_clump_mask.py')

if 'ind_map' not in dir():
    ind_map = ((na.mgrid[0:twod_clump.shape[0], 0:twod_clump.shape[0]] + 0.5)/twod_clump.shape[0]) #or some shit

clump_ind = 1
these_x = ind_map[0][twod_clump == clump_ind].flatten()
these_y = ind_map[1][twod_clump == clump_ind].flatten()
Left = na.zeros(3)
Right = na.ones(3)
nx, ny = twod_clump.shape
Left[ x_dict[0] ]  = na.floor(these_x.min()*nx)/nx
Left[ y_dict[0] ]  = na.floor(these_y.min()*nx)/nx
Right[ x_dict[0] ] = na.ceil(these_x.max()*ny)/ny
Right[ y_dict[0] ] = na.ceil(these_y.max()*ny)/ny 
Center = 0.5*(Left+Right)                         
region = ds.region(Center,Left,Right)
region.set_field_parameter('clump_mask_stuff',clump_stuff_tuple)
cut_region = region.cut_region(['obj["clump_mask"].astype("int") == %d'%clump_ind])
cut_region.set_field_parameter('clump_mask_stuff',clump_stuff_dict['clump_mask_stuff'])
uniqua =  na.unique(cut_region['clump_mask']) 
if 0:
    print uniqua
    if len(uniqua) != 1:
        print "ERROR with cut region, clump_ind = %d"%clump_ind,uniqua
