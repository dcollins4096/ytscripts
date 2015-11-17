
n_rows = max([len(gpd.name),len(gpc.name), len(gpf.name)])
n_column = max([len(n) for n in gpd.name+gpc.name+gpf.name])
n_column += 4
template_1 = "%"+str(n_column)+"s, "

for n in range(n_rows):
    for doer in [gpd,gpc,gpf]:
        mlist = doer.name
        if len(mlist) > n:
            print template_1%mlist[n],"%0.2f"%doer.self[n],
        else:
            print template_1%"--","----",
    print ""

