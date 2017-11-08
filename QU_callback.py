import yt
from yt.visualization.plot_modifications import *
import numpy as np
class PolFracCallback(PlotCallback):
    """
    Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints
    *scale* is the data units per arrow length unit using *scale_units*
    (see matplotlib.axes.Axes.quiver for more info)
    """
    _type_name = "polarized_angles"
    _supported_geometries = ("cartesian", "spectral_cube")
    def __init__(self, length = 1.0, factor=16, scale=None, weight=None,
                 scale_units=None, normalize=False, n0=1, p=1):
        PlotCallback.__init__(self)
        self.factor = factor
        self.scale = scale
        self.scale_units = scale_units
        self.normalize = normalize
        self.length = length
        self.weight=weight
        self.n0=n0
        self.p=p

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        # See the note about rows/columns in the pixelizer for more information
        # on why we choose the bounds we do
        nx = plot.image._A.shape[1] // self.factor
        ny = plot.image._A.shape[0] // self.factor
        # periodicity
        ax = plot.data.axis
        ds = plot.data.ds
        (xi, yi) = (ds.coordinates.x_axis[ax],
                    ds.coordinates.y_axis[ax])
        period_x = ds.domain_width[xi]
        period_y = ds.domain_width[yi]
        periodic = int(any(ds.periodicity))
        axis_dict={0:'x',1:'y',2:'z','x':'x','y':'y','z':'z'}
        Uname = "U%s_n0-%04d_p-%d"%(axis_dict[ax], self.n0, self.p)
        Qname = "Q%s_n0-%04d_p-%d"%(axis_dict[ax], self.n0, self.p)
        N2name = "N2%s_n0-%04d_p-%d"%(axis_dict[ax], self.n0, self.p)
        print "Uname = ", Uname
        print "Qname = ", Qname
        print "N2name = ", N2name
        #length = np.ones_like(theta_stokes)*self.length
        Q = plot.data[Qname]
        U = plot.data[Uname]
        N = plot.data['density'].in_units('code_density').v
        N2 = plot.data[N2name]
        theta_stokes = 0.5*np.arctan2(U,Q)*N2
        p0 = 0.1

        frac = p0*np.sqrt(Q**2+U**2)/(N - p0*N2)
        print "FRAC", frac.min(), frac.max()
        print "Q  ", Q.min(), Q.max(), np.isnan(Q).sum()/(Q.size*1.0)
        print "U  ", U.min(), U.max()
        print "N  ", N.min(), N.max()
        print "N2 ", N2.min(), N2.max()

        fv_x = frac*np.cos(theta_stokes)
        fv_y = frac*np.sin(theta_stokes)
        pixX = np.zeros((ny, nx), dtype="f8")
        pixY = np.zeros((ny, nx), dtype="f8")
        pixelize_cartesian(pixX, plot.data['px'], plot.data['py'],
                                  plot.data['pdx'], plot.data['pdy'],
                                  fv_x,
                                  (x0, x1, y0, y1), 0, # bounds, antialias
                                  (period_x, period_y), periodic)
        pixelize_cartesian(pixY, plot.data['px'], plot.data['py'],
                                  plot.data['pdx'], plot.data['pdy'],
                                  fv_y,
                                  (x0, x1, y0, y1), 0, # bounds, antialias
                                  (period_x, period_y), periodic)
        X,Y = np.meshgrid(np.linspace(xx0,xx1,nx,endpoint=True),
                          np.linspace(yy0,yy1,ny,endpoint=True))
        if self.normalize:
            nn = np.sqrt(pixX**2 + pixY**2)
            pixX /= nn
            pixY /= nn
        plot._axes.quiver(X,Y, pixX, pixY, scale=self.scale, scale_units=self.scale_units, headlength=0, headwidth=1)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
class MagdirCallback(PlotCallback):
    """
    Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints
    *scale* is the data units per arrow length unit using *scale_units*
    (see matplotlib.axes.Axes.quiver for more info)
    """
    _type_name = "magnetic_angles"
    _supported_geometries = ("cartesian", "spectral_cube")
    def __init__(self, length = 1.0, factor=16, scale=None,
                 scale_units=None, normalize=False, n0=1, p=1):
        PlotCallback.__init__(self)
        self.factor = factor
        self.scale = scale
        self.scale_units = scale_units
        self.normalize = normalize
        self.length = length
        self.n0=n0
        self.p=p

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        # See the note about rows/columns in the pixelizer for more information
        # on why we choose the bounds we do
        nx = plot.image._A.shape[1] // self.factor
        ny = plot.image._A.shape[0] // self.factor
        # periodicity
        ax = plot.data.axis
        ds = plot.data.ds
        (xi, yi) = (ds.coordinates.x_axis[ax],
                    ds.coordinates.y_axis[ax])
        period_x = ds.domain_width[xi]
        period_y = ds.domain_width[yi]
        periodic = int(any(ds.periodicity))
        (xi, yi) = (plot.data.ds.coordinates.x_axis[ax],
                    plot.data.ds.coordinates.y_axis[ax])
        axis_names = plot.data.ds.coordinates.axis_name
        xv = "magnetic_field_%s" % (axis_names[xi])
        yv = "magnetic_field_%s" % (axis_names[yi])
        B1 = plot.data[xv]
        B2 = plot.data[yv]
        this_norm = 1./np.sqrt(B1**2+B2**2)

        fv_x = this_norm*B1#self.length*np.cos(theta_stokes)
        fv_y = this_norm*B2#self.length*np.sin(theta_stokes)
        pixX = np.zeros((ny, nx), dtype="f8")
        pixY = np.zeros((ny, nx), dtype="f8")
        pixelize_cartesian(pixX, plot.data['px'], plot.data['py'],
                                  plot.data['pdx'], plot.data['pdy'],
                                  fv_x,
                                  (x0, x1, y0, y1), 0, # bounds, antialias
                                  (period_x, period_y), periodic)
        pixelize_cartesian(pixY, plot.data['px'], plot.data['py'],
                                  plot.data['pdx'], plot.data['pdy'],
                                  fv_y,
                                  (x0, x1, y0, y1), 0, # bounds, antialias
                                  (period_x, period_y), periodic)
        X,Y = np.meshgrid(np.linspace(xx0,xx1,nx,endpoint=True),
                          np.linspace(yy0,yy1,ny,endpoint=True))
        if self.normalize:
            nn = np.sqrt(pixX**2 + pixY**2)
            pixX /= nn
            pixY /= nn
        plot._axes.quiver(X,Y, pixX, pixY, scale=self.scale, scale_units=self.scale_units, headlength=0, headwidth=1, color='b')
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)

class StokesCallback(PlotCallback):
    """
    Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints
    *scale* is the data units per arrow length unit using *scale_units*
    (see matplotlib.axes.Axes.quiver for more info)
    """
    _type_name = "stokes_angles2"
    _supported_geometries = ("cartesian", "spectral_cube")
    def __init__(self, length = 1.0, factor=16, scale=None, weight=None,
                 scale_units=None, normalize=False, n0=1, p=1, kludge_my_package=None):
        PlotCallback.__init__(self)
        self.factor = factor
        self.scale = scale
        self.scale_units = scale_units
        self.normalize = normalize
        self.length = length
        self.weight=weight
        self.n0=n0
        self.p=p
        self.kludge_my_package = kludge_my_package

    def __call__(self, plot):
        print "call 1"
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        # See the note about rows/columns in the pixelizer for more information
        # on why we choose the bounds we do
        nx = plot.image._A.shape[1] // self.factor
        ny = plot.image._A.shape[0] // self.factor
        # periodicity
        ax = plot.data.axis
        ds = plot.data.ds
        (xi, yi) = (ds.coordinates.x_axis[ax],
                    ds.coordinates.y_axis[ax])
        period_x = ds.domain_width[xi]
        period_y = ds.domain_width[yi]
        periodic = int(any(ds.periodicity))
        axis_dict={0:'x',1:'y',2:'z','x':'x','y':'y','z':'z'}
        Uname = "U%s_n0-%04d_p-%d"%(axis_dict[ax], self.n0, self.p)
        Qname = "Q%s_n0-%04d_p-%d"%(axis_dict[ax], self.n0, self.p)
        print "Uname = ", Uname
        print "Qname = ", Qname
        #length = np.ones_like(theta_stokes)*self.length
        Q = plot.data[Qname]
        U = plot.data[Uname]
        I = np.sqrt(Q**2+U**2)
        theta_stokes = 0.5*np.arctan2(U,Q)
        norm = self.length
        if self.weight == 'N':
            N = plot.data['density']
            norm = N/N.max()
        elif self.weight=="I":
            norm = I #/I.max()
            self.scale = self.scale/I.max().v
            print "Weight I (rly)"
        print "WEIGHT WTF hax2", self.weight
        #theta_stokes = np.pi/7.
        fv_x = norm*np.cos(theta_stokes)
        fv_y = norm*np.sin(theta_stokes)
        pixX = np.zeros((ny, nx), dtype="f8")
        pixY = np.zeros((ny, nx), dtype="f8")
        pixelize_cartesian(pixX, plot.data['px'], plot.data['py'],
                                  plot.data['pdx'], plot.data['pdy'],
                                  fv_x,
                                  (x0, x1, y0, y1), 0, # bounds, antialias
                                  (period_x, period_y), periodic)
        pixelize_cartesian(pixY, plot.data['px'], plot.data['py'],
                                  plot.data['pdx'], plot.data['pdy'],
                                  fv_y,
                                  (x0, x1, y0, y1), 0, # bounds, antialias
                                  (period_x, period_y), periodic)
        X,Y = np.meshgrid(np.linspace(xx0,xx1,nx,endpoint=True),
                          np.linspace(yy0,yy1,ny,endpoint=True))
        #global kludge_my_package
        self.kludge_my_package['Q']=Q
        self.kludge_my_package['U']=U
        self.kludge_my_package['theta_stokes']=theta_stokes
        self.kludge_my_package['X']=X
        self.kludge_my_package['Y']=Y
        self.kludge_my_package['pixX']=pixX
        self.kludge_my_package['pixY']=pixY
        self.kludge_my_package['fv_x']=fv_x
        self.kludge_my_package['fv_y']=fv_y
        self.kludge_my_package['I']=I
        self.kludge_my_package['norm']=norm
        if self.normalize:
            nn = np.sqrt(pixX**2 + pixY**2)
            pixX /= nn
            pixY /= nn
        plot._axes.quiver(X,Y, pixX, pixY, scale=self.scale, scale_units=self.scale_units, headlength=0, headwidth=1,
                          color='r')
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        print "YES I RELOADED"


class StokesVectorsCallback_a(PlotCallback):
    """
    Adds a 'quiver' plot of velocity to the plot, skipping all but
    every *factor* datapoint. *scale* is the data units per arrow
    length unit using *scale_units* (see
    matplotlib.axes.Axes.quiver for more info). if *normalize* is
    True, the velocity fields will be scaled by their local
    (in-plane) length, allowing morphological features to be more
    clearly seen for fields with substantial variation in field
    strength.
    """
    _type_name = "stokes_vectors_take1"
    _supported_geometries = ("cartesian", "spectral_cube")
    def __init__(self, factor=16, scale=None, scale_units=None, normalize=False, n0=1, p=1, length=0.5):
        PlotCallback.__init__(self)
        self.factor = factor
        self.scale  = scale
        self.scale_units = scale_units
        self.normalize = normalize
        self.length=length
        self.n0 = n0
        self.p=p

    def __call__(self, plot):
        # Instantiation of these is cheap
        if plot._type_name == "CuttingPlane":
            raise
            qcb = CuttingQuiverCallback("cutting_plane_velocity_x",
                                        "cutting_plane_velocity_y",
                                        self.factor, scale=self.scale,
                                        normalize=self.normalize,
                                        scale_units=self.scale_units)
        else:
            ax = plot.data.axis
            (xi, yi) = (plot.data.ds.coordinates.x_axis[ax],
                        plot.data.ds.coordinates.y_axis[ax])
            axis_names = plot.data.ds.coordinates.axis_name
            U = "U%s_n0-%04d_p-%d"%(axis_names[x], self.n0, self.p)
            Q = "Q%s_n0-%04d_p-%d"%(axis_names[x], self.n0, self.p)
            theta_stokes = 0.5*np.arctan2(U,Q)
            length = np.ones_like(theta_stokes)*self.length
            xv = length*np.cos(theta)
            yv = length*np.sin(theta)

            qcb = QuiverCallback(xv, yv, self.factor, scale=self.scale,
                                 scale_units=self.scale_units,
                                 normalize=self.normalize, bv_x=bv_x, bv_y=bv_y)
        return qcb(plot)

class ThetaQuiverCallback(PlotCallback):
    """
    Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints
    *scale* is the data units per arrow length unit using *scale_units*
    (see matplotlib.axes.Axes.quiver for more info)
    """
    _type_name = "theta_quiver"
    _supported_geometries = ("cartesian", "spectral_cube")
    def __init__(self, theta, length, factor=16, scale=None,
                 scale_units=None, normalize=False):
        PlotCallback.__init__(self)
        self.theta = theta
        self.length = length
        self.factor = factor
        self.scale = scale
        self.scale_units = scale_units
        self.normalize = normalize

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        # See the note about rows/columns in the pixelizer for more information
        # on why we choose the bounds we do
        nx = plot.image._A.shape[1] // self.factor
        ny = plot.image._A.shape[0] // self.factor
        # periodicity
        ax = plot.data.axis
        ds = plot.data.ds
        (xi, yi) = (ds.coordinates.x_axis[ax],
                    ds.coordinates.y_axis[ax])
        period_x = ds.domain_width[xi]
        period_y = ds.domain_width[yi]
        periodic = int(any(ds.periodicity))
        fv_x = plot.data[self.field_x]
        if self.bv_x != 0.0:
            # Workaround for 0.0 without units
            fv_x -= self.bv_x
        fv_y = plot.data[self.field_y]
        if self.bv_y != 0.0:
            # Workaround for 0.0 without units
            fv_y -= self.bv_y
        pixX = np.zeros((ny, nx), dtype="f8")
        pixY = np.zeros((ny, nx), dtype="f8")
        pixelize_cartesian(pixX, plot.data['px'], plot.data['py'],
                                  plot.data['pdx'], plot.data['pdy'],
                                  fv_x,
                                  (x0, x1, y0, y1), 0, # bounds, antialias
                                  (period_x, period_y), periodic)
        pixelize_cartesian(pixY, plot.data['px'], plot.data['py'],
                                  plot.data['pdx'], plot.data['pdy'],
                                  fv_y,
                                  (x0, x1, y0, y1), 0, # bounds, antialias
                                  (period_x, period_y), periodic)
        X,Y = np.meshgrid(np.linspace(xx0,xx1,nx,endpoint=True),
                          np.linspace(yy0,yy1,ny,endpoint=True))
        if self.normalize:
            nn = np.sqrt(pixX**2 + pixY**2)
            pixX /= nn
            pixY /= nn
        plot._axes.quiver(X,Y, pixX, pixY, scale=self.scale, scale_units=self.scale_units)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
