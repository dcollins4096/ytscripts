
proj = ds.proj('density',0,data_source=leaf_clumps[0].data,center=position)
pw = proj.to_pw(center=position,width = 0.1)
print pw.save('p14b_atest_leaf')
