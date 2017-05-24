if 'ef' not in dir():
    execfile('go')

flt = taxi.fleet(['g01','g02'])
#flt['fields'] = ['Metal_Density'] + MultiSpecies1
flt['fields']=['density']
flt['callbacks']=['star_particles']
flt['particle_filter_names']=['formed_star']
flt.plot(prefix='p33_')

