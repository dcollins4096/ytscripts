execfile('go')
execfile('p42_pg.py')
execfile('p42_spectra.py')
ef('p42_quan.py')
#aw14=taxi.taxi(directory='/scratch/00369/tg456484/Paper49_EB/aw14_M3_MA0.3_256_p59',name='aw14',frames=[100,150])
#aw12=taxi.taxi(directory='/scratch/00369/tg456484/Paper49_EB/aw12_M1_Ma0.3_256_p59',name='aw12',frames=[100,150])
#aw13=taxi.taxi(directory='/scratch/00369/tg456484/Paper49_EB/aw13_M3_MA1_256_p59',name='aw13',frames=[100,150,200])
aq32 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq32_m2.9_drive1_32',name='aq32',frames=range(0,42,2),fields=['density'])
aq31 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq31_m2.9_drive0.5_32',name='aq31',frames=range(0,42,2),fields=['density'])
aq16 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq16_ppm_m9_drive0_noamr_128',name='aq15',frames=range(17),fields=['density'])
aq33 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq33_m9_drive0.5_32',name='aq33',frames=range(0,110,10))
aq34 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq34_m9_drive0_p60_32', name='aq34', frames=range(0,110,10))
aq38 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq38_m2.9_drive0_32_p59_ppm',name='aq38',frames=range(0,42,2),fields=['density'])
aq39 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq39_m2.9_drive0_32_p59_ppm' ,name='aq39',frames=range(0,42,2),fields=['density'])
aq40 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq40_m2.9_drive0_32_p59_hydro6_newAddVel' ,name='aq40',frames=range(0,42,2),fields=['density'])
aq41 = taxi.taxi(directory='/scratch1/dcollins/Paper42_new_turb/aq41_m9_drive0.5_32_p59_ppm' ,name='aq41',frames=range(0,42,2),fields=['density'])
aq42 = taxi.taxi(directory='/scratch/00369/tg456484/Paper42_NewAK/aq42_m9_drive0.5_512_p59_ppm' ,name='aq42',frames=range(90,0,-10),fields=['density'])
aw11 = taxi.taxi(directory='/scratch/00369/tg456484/Paper49_EB/aw11_M1_MA1_256_p59' ,name='aw11',frames=range(0,210,10),fields=['density'])
#aq44 = taxi.taxi(directory='/scratch/00369/tg456484/Paper42_NewAK/aq44_Actually9_d0.5_p59',name='aq44',frames=[int(sys.argv[1])])
aw15 = taxi.taxi('aw15')
aw16 = taxi.taxi('aw16')
aw17 = taxi.taxi('aw17')
aw18 = taxi.taxi('aw18')
fields = ['%s-velocity'%s for s in 'xyz']
#fields = ['DrivingField%s'%s for s in '123']
#all_the_spectra(aw14,fields)
#all_the_spectra(aw12,fields)
#all_the_spectra(aw13,fields)
#all_the_spectra(aq44,fields)
fleet = [aw15,aw16,aw17,aw18]
quan_boxes = {}
for car in fleet:
    car.fill(0)
    all_frames = car.frame_dict.keys()
    car.frames = all_frames[::10]
    if all_frames[-1] not in car.frames:
        car.frames += [all_frames[-1]]
    quan_boxes[car.name] = quan_box(car)
    quan_boxes[car.name](car)
#for car in fleet:
#car.plot()
