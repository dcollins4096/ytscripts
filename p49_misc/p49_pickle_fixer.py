import fPickle
sim='b20_512'
a = fPickle.load('quan_box_%s_cmbtn.pickle'%sim)
b = fPickle.load('quan_box_%s.pickle'%sim)
new_thing = b
new_thing.pop('EB')
new_thing['EB']=a['EB']
fPickle.dump(new_thing,'quan_box_%s_both.pickle'%sim)

