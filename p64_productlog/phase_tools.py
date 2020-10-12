from go import *
class phase_things():

    def __init__(self,phase,field='cell_volume'):
        self.phase=phase
        self.field=field

    #
    # Set up the coordinate fields.
    #
        x_bins = phase.x_bins
        y_bins = phase.y_bins
        x_cen = 0.5*(x_bins[1:]+x_bins[:-1])
        y_cen = 0.5*(y_bins[1:]+y_bins[:-1])
        nbins1 = x_cen.size
        nbins2 = y_cen.size
        #this line turns 1d coordinate arrays x_cen & ycen into  2d coordinate fields.
        self.TheX = np.r_[(nbins2)*[x_cen]].transpose()
        self.TheY = np.r_[(nbins1)*[y_cen]]

        self.x_del = (x_bins[1:]-x_bins[:-1])
        self.y_del = (y_bins[1:]-y_bins[:-1])
        self.x_del.shape = (self.x_del.size,1)
        self.dv = 1./(self.x_del*self.y_del)

        #The quantity to plot
        self.Pjoint = phase[field]
        self.CCC = self.Pjoint*self.dv

        self.P_x =(self.Pjoint).sum(axis=1)
        self.P_y =(self.Pjoint).sum(axis=0)
        self.P_x.shape = (self.P_x.size,1)
        self.IndJoint = self.P_x*self.P_y*self.dv

        #self.P_x =(self.Pjoint).sum(axis=1)
        #self.P_y =(self.Pjoint).sum(axis=0)
        #self.P_x.shape = (self.P_x.size,1)
        #IndJoint = self.P_x*self.P_y*self.dv

    def image(self,fig,ax):
        #Matplotlib color wizardry.
        vmin = self.CCC[self.CCC>0].min()
        vmax = self.CCC.max()
        print(vmin,vmax)
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
        cmap = mpl.cm.jet
        cmap.set_under('w')
        color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

        #Color field plot
        ploot=ax.pcolormesh(self.TheX,self.TheY,self.CCC,norm=norm)
        ax.set_xlabel(self.phase.x_field)
        ax.set_ylabel(self.phase.y_field)
        fig.colorbar(ploot,ax=ax)

    def image_joint(self,fig,ax):
        #
        # Also compare the outer product of the marginalized
        # distributions.  Take 1, skipping the bin size.
        #
        #pdb.set_trace()
        TheZ = self.IndJoint
        vmin = TheZ[TheZ>0].min()
        vmax = TheZ.max()
        print(vmin,vmax)
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
        plot2=ax.pcolormesh(self.TheX,self.TheY,TheZ,norm=norm)
        ax.set_xlabel(self.phase.x_field)
        ax.set_ylabel(self.phase.y_field)
        fig.colorbar(plot2,ax=ax)

    #
    # Lets also try it with contours.
    #
    def plot_contours(self,fig,ax):
        levels=np.arange(0,self.CCC.max(),self.CCC.max()/7)
        ax.contour(self.TheX,self.TheY,self.CCC,levels,colors='k')
        ax.contour(self.TheX,self.TheY,self.IndJoint,levels,colors='r')
        ax.set_xlabel(self.phase.x_field)
        ax.set_ylabel(self.phase.y_field)


def plot_phase_contours_1(fig,ax,phase, field='cell_volume', ax2=None):


    #
    # Set up the coordinate fields.
    #
    x_bins = phase.x_bins
    y_bins = phase.y_bins
    x_cen = 0.5*(x_bins[1:]+x_bins[:-1])
    y_cen = 0.5*(y_bins[1:]+y_bins[:-1])
    nbins1 = x_cen.size
    nbins2 = y_cen.size
    #this line turns 1d coordinate arrays x_cen & ycen into  2d coordinate fields.
    TheX = np.r_[(nbins2)*[x_cen]].transpose()
    TheY = np.r_[(nbins1)*[y_cen]]

    x_del = (x_bins[1:]-x_bins[:-1])
    y_del = (y_bins[1:]-y_bins[:-1])
    x_del.shape = (x_del.size,1)
    dv = 1./(x_del*y_del)

    #The quantity to plot
    Pjoint = phase[field]
    CCC = Pjoint*dv

    #Matplotlib color wizardry.
    norm = mpl.colors.Normalize(vmin=CCC[CCC>0].min(),vmax=CCC.max())
    cmap = mpl.cm.jet
    cmap.set_under('w')
    color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

    #Color field plot
    ploot=ax.pcolormesh(TheX,TheY,CCC,norm=norm)
    ax.set_xlabel(phase.x_field)
    ax.set_ylabel(phase.y_field)
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.colorbar(ploot,ax=ax)

    if ax2:
        #
        # Also compare the outer product of the marginalized
        # distributions.  Take 1, skipping the bin size.
        #
        P_x =(Pjoint).sum(axis=1)
        P_y =(Pjoint).sum(axis=0)
        P_x.shape = (P_x.size,1)
        IndJoint = P_x*P_y*dv
        #pdb.set_trace()
        plot2=ax2.pcolormesh(TheX,TheY,IndJoint,norm=norm)
        ax2.set_xlabel(phase.x_field)
        ax2.set_ylabel(phase.y_field)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        fig.colorbar(plot2,ax=ax2)

    #
    # Lets also try it with contours.
    #
    fig6, ax6 = plt.subplots(1,1)
    levels=np.arange(0,CCC.max(),CCC.max()/7)
    ax6.contour(TheX,TheY,CCC,levels,colors='k')
    ax6.contour(TheX,TheY,IndJoint,levels,colors='r')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_xlabel(phase.x_field)
    ax6.set_ylabel(phase.y_field)
    fig6.savefig('../PigPen/contour_test.png')
    plt.close(fig6)
