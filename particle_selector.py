dbg = 1
import pyximport; pyximport.install()
import particle_ops

def deposit_target_particles_1(field, data):
    indices_late = data.get_field_parameter('indices_late')
    blank = np.zeros(data.ActiveDimensions, dtype='float64')
    if data.NumberOfParticles == 0: return blank
    if dbg > 0:
        if hasattr(indices_late,'sum'):
            print "indices inside", indices_late.sum()
        else:
            print "indices inside", indices_late
            return blank
    #int 32 because it's just a flag.
    #mask_to_get = np.zeros(indices_late.shape, dtype='int32')
    mask_to_get = data.get_field_parameter('mask_to_get')
    my_indices = data['particle_index'].astype('int64')
    print "  min thing"
    #mask_to_get[ indices_late < my_indices.min()] = 1
    #mask_to_get[ indices_late > my_indices.max()] = 1
    print "  left", mask_to_get.size - mask_to_get.sum()
    print "  search"
    found_any, mask = particle_ops.mask_particles(
        indices_late, my_indices, mask_to_get)
    data.set_field_parameter('mask_to_get',mask_to_get)
    pos_to_get = data['particle_position'][mask == 1]
    d = data.deposit(pos_to_get, method = "count")
    d = data.ds.arr(d, input_units = "cm**-3")
    return data.apply_units(d, field.units)

yt.add_field(      ("deposit","deposit_target_particles_1"),
         function = deposit_target_particles_1,
         validators = [yt.ValidateSpatial(), 
                       yt.ValidateParameter('indices_late'), 
                       yt.ValidateParameter('mask_to_get'), 
                       yt.ValidateGridType()],
         display_name = "target_particles")

