from go import *
import turb_quan
reload(turb_quan)
reload(taxi)
ab22=taxi.taxi('ab22')
qb = turb_quan.quan_box(ab22)

qb.load(h5_name = "p49_data_h5/quan_box_%s.h5"%'ab22')
