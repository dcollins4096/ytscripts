from go import *

ddd = "/scratch/00369/tg456484/g47_xsede_2022/Agora_Full"
frame=100
sim = "%s/DD%04d/DD%04d"%(ddd,frame,frame)
ds = yt.load(sim)
