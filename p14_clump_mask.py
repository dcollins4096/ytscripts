#execfile('p14_start.py')

def region_projection(field,data):
    #global flat_clump, clump_axis,twod_clump
    flag_array = np.zeros_like(data['density']).v -1
    if data.has_field_parameter('clump_mask_stuff'):
        clump_stuff_tuple = data.get_field_parameter('clump_mask_stuff')
        if type(clump_stuff_tuple) is tuple:
            (flat_clump,clump_axis,twod_clump)=clump_stuff_tuple 
            x_dict = [1,0,0]
            y_dict = [2,2,1]
            x_like = 'xyz'[x_dict[clump_axis]]
            y_like = 'xyz'[y_dict[clump_axis]]
            clump_dx = 1./twod_clump.shape[0]
            clump_dy = 1./twod_clump.shape[1]
            index_x = (data[x_like]/clump_dx).astype('int').flatten()
            index_y = (data[y_like]/clump_dy).astype('int').flatten()
            #index_array = (data[y_like]/clump_dy + clump.shape[1]*data[x_like]/clump_dx).astype('int')
            index_array = (index_y + twod_clump.shape[1]*index_x).astype('int')
            flag_array = flat_clump[index_array].astype(data[x_like].dtype)
            flag_array.shape = data[x_like].shape
    return flag_array
#clump_mask_defined = True
fieldname='clump_mask'
yt.add_field(fieldname,function=region_projection,validators=[yt.ValidateSpatial(0)],take_log=False)

#oober.h.grids[0][fieldname]
#oober.fields=[fieldname]
#oober.plot()
#img_region(region,fields=['region_projection'],axis=[0,1,2],prefix='p14_region_test')
