

def plot_clump(clump, ax=[0], field='density',cmap=None, contour=True, bump=False, plt_args = {}, prefix = 'porkchop'):
    #get the ds from the parent chain.
    tmp = clump
    while tmp.parent is not None:
        tmp = tmp.parent
    ds = tmp.data.ds
    x = clump['x']; y = clump['y']; z = clump['z']
    Left  = nar([x.min(), y.min(), z.min()])
    Right = nar([x.max(), y.max(), z.max()])
    Center = 0.5*(Left + Right)
    if bump:
        for i in [0,1,2]:
            Left[i] = max(min(Left[i], Center[i]-0.05), 0)
            Right[i] = min(max(Right[i], Center[i]+0.05), 1)


    Width = Right - Left
    region = ds.region(Center, Left, Right)
    for axis in ax:
        PlotWidth = max(Width[axis-1], Width[axis-2]) #I'm quite pleased with this.
        proj = ds.proj( field, axis, data_source = region, center = Center)
        pw = proj.to_pw(center=Center,width=(PlotWidth, 'code_length'))
        pw.annotate_clumps([clump])
        pw.save(prefix)

for n, cl in enumerate(leaf_clumps[n:]):
    outname = 'u05_0125_cl%04d'%n
    plot_clump(cl,prefix=outname, bump=True)
