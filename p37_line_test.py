

dx=1.; dy=1.; NpointsX = 20.*1j; NpointsY=30.*1j
rightx = 20.; righty = 30.
y, x = np.mgrid[0.5*dx:rightx-0.5*dx:NpointsX, 0.5*dy:righty-0.5*dy:NpointsY]


def point_selector(x,y,x0,x1,y0,y1,dx,dy):
    #F(x,y) = (y1-y0)*x + (x0-x1)*y +(x1*y0-x0*y1)
    #F==0, (x,y) is on the line.  <0 above, >0 below
    #All above or below, then there's no intersection. 
    #Sign change means that two points are on opposite sides of the line.

    c1 = y1-y0
    c2 = x0-x1
    c3 = x1*y0-x0*y1
    a   = c1*(x      )+c2*(y      )+c3
    amm = c1*(x-dx/2.)+c2*(y-dy/2.)+c3
    amp = c1*(x-dx/2.)+c2*(y+dy/2.)+c3
    app = c1*(x+dx/2.)+c2*(y+dy/2.)+c3
    apm = c1*(x+dx/2.)+c2*(y-dy/2.)+c3
    keep = amm*amp <= 0
    keep = np.logical_or( keep, amp*app <= 0 )
    keep = np.logical_or( keep, app*apm <= 0 )
    keep = np.logical_or( keep, apm*amm <= 0 ) 
    #pdb.set_trace()
    return keep

xf0=5.; xf1 = 7.
yf0=0.; yf1 = 20.0
keep= point_selector(x,y,xf0,xf1,yf0,yf1,dx,dy)
print keep.astype('int')
def prnt(arr,frm="  %5.1f"):
    nx,ny = arr.shape
    for i in range(ny):
        print frm*nx%tuple(arr[i])
#prnt(keep)

"""
y = m x + b
m= (y1-y0)/(x1-x0)
 (x1-x0) y - (y1-y0) x =
 (x1-x0)y1 - (y1-y0)x1 = b
 x1 y1-x0 y1 - y1 x1 + y0 x1 = b
 ----          -----
 y0x1 - x0y1 = b
 0 = -(x1-x0)y + (y1-y0)x + b
 """
