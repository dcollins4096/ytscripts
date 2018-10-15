execfile('go')
base = '/scratch/00369/tg456484/Paper49/'
name = 'aa19_M1_Ma0.3_256'
fields=['density','magnetic_field_strength']
frames='every 10'

aa19 = taxi.taxi(dir=base+name,name='aa19',fields=fields,frames=frames)
aa20 = taxi.taxi(dir=base+'aa20_M3_Ma1_256', name = 'aa20',fields=fields,frames=frames)
aa22 = taxi.taxi(dir=base+'aa21_M0.6_MA0.3_256', name = 'aa21',fields=fields,frames=frames)
aa22 = taxi.taxi(dir=base+'aa22_M3_Ma03_256', name = 'aa22',fields=fields,frames=frames)
az19 = taxi.taxi(dir=base+'az19_M1_Ma0.3_512_cfl', name = 'az19',fields=fields,frames=frames)
az20 = taxi.taxi(dir=base+'az20_M3_Ma1_512', name = 'az20',fields=fields,frames=frames)
az21 = taxi.taxi(dir=base+'az21_M0.6_MA0.3_512', name = 'az21',fields=fields,frames=frames)
az22 = taxi.taxi(dir=base+'az22_M3_Ma03_512', name = 'az22',fields=fields,frames=frames)

aa19.save('aa19')
aa20.save('aa20')
aa22.save('aa22')
aa22.save('aa22')
az19.save('az19')
az20.save('az20')
az21.save('az21')
az22.save('az22')
