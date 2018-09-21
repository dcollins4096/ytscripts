#from yt.utilities.lib import CICDeposit_3
def _particle_image_1(field, data):
    global ybucket, zbucket
    """Returns 1 where particles are, 0 where they are not."""
    blank = na.zeros(data.ActiveDimensions, dtype='float64')
    if data.NumberOfParticles == 0: return blank
    indices_late = data.get_field_parameter('indices_late')
    print indices_late.sum()
    #mask applied to indices_late: stops redudnant searching
    if 0:
        mask_to_get = na.zeros(indices_late.shape, dtype='int32')
        found_any, mask = particle_ops.mask_particles(
            indices_late, data['particle_index'], mask_to_get)
        xpos = data['particle_position_x'][ mask == 1]
        ypos = data['particle_position_y'][ mask == 1]
        zpos = data['particle_position_z'][ mask == 1]
        ybucket.append(ypos)
        zbucket.append(zpos)
        NMaskedParticles = na.int64(xpos.size)
        unit_mass = na.ones(NMaskedParticles,dtype = "float64")
        CICDeposit_3(xpos.astype(na.float64),
                     ypos.astype(na.float64),
                     zpos.astype(na.float64),
                     unit_mass,
                     NMaskedParticles,
                     blank, 
                     na.array(data.LeftEdge).astype(na.float64),
                     na.array(data.ActiveDimensions).astype(na.int32),
                     na.float64(data['dx']))
    min_density = 0.5*data.dds.prod()
    #blank[ blank > min_density ] = 1
    return blank
validator_list =[yt.ValidateParameter('indices_late'),yt.ValidateGridType()]
yt.add_field(         "particle_image_1", 
          function=_particle_image_1,
          validators=validator_list,take_log=False)
