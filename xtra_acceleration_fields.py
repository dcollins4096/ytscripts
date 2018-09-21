if 0:
    sl_center = slice(1, -1, None)
    sl_left = slice(None, -2, None)
    sl_right = slice(2, None, None)
    div_fac = 2.0
    ftype = 'gas'
    def _Avorticity_x(field, data):
        f  = (data[ftype, "acceleration_z"][sl_center,sl_right,sl_center] -
              data[ftype, "acceleration_z"][sl_center,sl_left,sl_center]) \
              / (div_fac*yt.funcs.just_one(data["index", "dy"]).in_cgs())
        f -= (data[ftype, "acceleration_y"][sl_center,sl_center,sl_right] -
              data[ftype, "acceleration_y"][sl_center,sl_center,sl_left]) \
              / (div_fac*yt.funcs.just_one(data["index", "dz"].in_cgs()))
        new_field = data.ds.arr(np.zeros_like(data[ftype, "acceleration_z"],
                                              dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field
    def _Avorticity_y(field, data):
        f  = (data[ftype, "acceleration_x"][sl_center,sl_center,sl_right] -
              data[ftype, "acceleration_x"][sl_center,sl_center,sl_left]) \
              / (div_fac*yt.funcs.just_one(data["index", "dz"]))
        f -= (data[ftype, "acceleration_z"][sl_right,sl_center,sl_center] -
              data[ftype, "acceleration_z"][sl_left,sl_center,sl_center]) \
              / (div_fac*yt.funcs.just_one(data["index", "dx"]))
        new_field = data.ds.arr(np.zeros_like(data[ftype, "acceleration_z"],
                                              dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field
    def _Avorticity_z(field, data):
        f  = (data[ftype, "acceleration_y"][sl_right,sl_center,sl_center] -
              data[ftype, "acceleration_y"][sl_left,sl_center,sl_center]) \
              / (div_fac*yt.funcs.just_one(data["index", "dx"]))
        f -= (data[ftype, "acceleration_x"][sl_center,sl_right,sl_center] -
              data[ftype, "acceleration_x"][sl_center,sl_left,sl_center]) \
              / (div_fac*yt.funcs.just_one(data["index", "dy"]))
        new_field = data.ds.arr(np.zeros_like(data[ftype, "acceleration_z"],
                                              dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field

    Avort_validators = [yt.ValidateSpatial(1,
                            [(ftype, "acceleration_x"),
                             (ftype, "acceleration_y"),
                             (ftype, "acceleration_z")])]
    name = []                    
    for ax in 'xyz':
        n = "Avorticity_%s" % ax

        yt.add_field((ftype, n),
                           function=eval("_%s" % n),
                           units="1/cm",
                           validators=Avort_validators)

    def create_magnitude_field(registry, basename, field_units,
                               ftype="gas", slice_info=None,
                               validators=None, particle_type=False):

        field_components = [(ftype, "%s_%s" % (basename, ax)) for ax in 'xyz']

        def _magnitude(field, data):
            fn = field_components[0]
            mag = data[fn] * data[fn]
            for idim in range(1, 3):
                fn = field_components[idim]
                mag += data[fn] * data[fn]
            return np.sqrt(mag)

        registry.add_field((ftype, "%s_magnitude" % basename),
                           function=_magnitude, units=field_units,
                           validators=validators, particle_type=particle_type)
    create_magnitude_field(yt, "Avorticity", "1/cm",
                           ftype=ftype, slice_info=None,
                           validators=Avort_validators)
    def _Adivergence(field, data):
        xn,yn,zn = ["acceleration_%s" % ax for ax in 'xyz']
        ds = div_fac * just_one(data["index", "dx"])
        f  = data[xn][sl_right,1:-1,1:-1]/ds
        f -= data[xn][sl_left ,1:-1,1:-1]/ds
        ds = div_fac * just_one(data["index", "dy"])
        f += data[yn][1:-1,sl_right,1:-1]/ds
        f -= data[yn][1:-1,sl_left ,1:-1]/ds
        ds = div_fac * just_one(data["index", "dz"])
        f += data[zn][1:-1,1:-1,sl_right]/ds
        f -= data[zn][1:-1,1:-1,sl_left ]/ds
        new_field = data.ds.arr(np.zeros(data[xn].shape, dtype=np.float64),
                                f.units)
        new_field[1:-1,1:-1,1:-1] = f
        return new_field

    def _divergence_abs(field, data):
        return np.abs(data[ftype, "%s_divergence" % basename])

    yt.add_field((ftype, "%s_divergence" % 'acceleration'),
                       function=_Adivergence,
                       units="1/cm",
                       validators=[yt.ValidateSpatial(1)])
