
#def _very_silly(field,data):
#    return data['Bx']
#yt.add_field('very_silly',  function=_very_silly, validators=[yt.ValidateGridType()], units='code_magnetic')

ds1=yt.load('/scratch1/dcollins/Paper19/SphereTest/s03_mhdct/DD0002/data0002')
print ds1.index.grids[-1]['very_silly']
ds2=yt.load('/scratch1/dcollins/Paper19/SphereTest/s02_full/DD0002/data0002')
print ds2.index.grids[-1]['density']

