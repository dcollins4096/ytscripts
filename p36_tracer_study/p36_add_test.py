
c01=taxi.taxi(directory='/scratch1/dcollins/Paper36_tracertests/AddPost/c01_small_sphere_particles',name='c01',frames=[0,1,20],fields=['density'],callbacks=['particles'])
c01.cmap['density']='gray'
#c01.plot()
c03=taxi.taxi(directory='/scratch1/dcollins/Paper36_tracertests/AddPost/c03_c02_added',name='c03',frames=[0,1,20],fields=['density'],callbacks=['particles'])
c03.cmap['density']='gray'
#c03.plot()
c04=taxi.taxi(directory='/scratch1/dcollins/Paper36_tracertests/AddPost/c04_sphere_amr',name='c04',frames=[2],fields=['density'],callbacks=['particles'])
c04.cmap['density']='gray'
c04.plot()
c06=taxi.taxi(directory='/scratch1/dcollins/Paper36_tracertests/AddPost/c06_c05_add',name='c06',frames=[0],fields=['density'],callbacks=['particles'])
c06.cmap['density']='gray'
c06.plot()
