from go import *
import turb_quan
reload(turb_quan)
reload(taxi)
car=taxi.load('bb_b1p.1')
qb = turb_quan.quan_box(car)
qb.load()#h5_name = "p49_data_h5/quan_box_%s.h5"%'ab22')
qb.GetQUEB(900)
qb.PlotQUEB(900)
