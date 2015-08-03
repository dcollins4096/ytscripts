
if 1:
    subset='phil_trace_1'
    L, R = np.zeros(3), np.zeros(3)
    off = [-0.025, 0.14]
    anchor={80:[0.5,1.5],70:[0.5,1.5],60:[0.5,1.0],50:[0.5,1.0]}
    trace={80:[-571.,-492.],70:[-632.,-584.],60:[-616.,0],50:[-616.,-149.]}
    L[0] = anchor[frame][0]/4.6+trace[frame][0]/8192 + off[0]
    L[1] = anchor[frame][1]/4.6 +trace[frame][1]/8192 + off[1]
    R[0] = L[0]+700./8192
    R[1] = L[1]+700./8192
    L[2]=0.0; R[2]=1.0
    C = 0.5*(L+R)
    print "L",L
    print "R",R
