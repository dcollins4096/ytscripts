if 'ef' not in dir():
    execfile('go')
if 'sphere' not in dir():
    dataset = '/scratch1/dcollins/Paper08/B02/256/RS0050/restart0050'
    ds = yt.load(dataset)
    #max_density, position = ds.find_max('Density')
    sphere = ds.sphere(position, 0.05)
    #cs = clump_stuff.clump_set(data_object=[sphere],c_min=10,c_max=sphere['Density'].max(),step_size=1.1,field='Density')

import clump_properties
reload(clump_properties)
b = clump_properties.clump_stuff(sphere,1)
