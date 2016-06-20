
if 0:
    def _scaled_div_b(field,data):
        sdb = np.abs(data['enzo','DivB'])
        sdb /= data['magnetic_field_strength']
        sdb *= data.dds.max()
        return sdb
    yt.add_field('scaled_div_b',  function=_scaled_div_b, validators=[yt.ValidateGridType()])

    def _abs_divb(field,data):
        return np.abs(data['DivB'])
    yt.add_field('abs_divb', function = _abs_divb)

#print "MONKEY ON THE CAR"
