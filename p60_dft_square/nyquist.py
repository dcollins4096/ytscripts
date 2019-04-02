
N1 = 100
x1 = np.arange(0,2*np.pi,2*np.pi/(10*N1))
x2 = np.arange(0,2*np.pi,2*np.pi/N1)
Ny = N1//2
output_dir = "p60_dft_square/plots/"
if 0:
    for Dny in range(1,5): #range(Ny-10,Ny+1,2):
        k=Ny-Dny
        plt.clf()
        y1 = np.cos(k*x1)
        y2 = np.cos(k*x2)
        if 0:
            #plot the originals showing the beat frequency
            plt.plot(x1,y1,c='r')
            plt.plot(x2,y2,c='g')
        if 0:
            #plot the weird ass result
            other_nyquist=np.sin((Dny)*x1)
            super_nyquist=0.5*(np.cos((Ny+Dny)*x1)+np.cos((Ny-Dny)*x1))
            plt.plot(x2,y2,c='g')
            plt.plot(x1,super_nyquist,c='k')
        if 0:
            #plot the approximation
            other_nyquist=np.sin((Dny)*x1)
            super_nyquist=0.5*(np.cos((Ny+Dny)*x1)+np.cos((Ny-Dny)*x1))
            plt.plot(x2,y2-super_nyquist[::10],c='k')
        #plt.plot(x1,other_nyquist,c='r')

        #plt.plot(x1,np.cos(Dny*x1),c='k')
        #plt.plot(x1,np.cos((Ny+Dny)*x1),c=[0.5]*3)
        #plt.plot(x1,,c=[0.5]*3)

        plt.savefig('%s/square_nyquist.pdf'%output_dir)
        time.sleep(1)


if 1:
    for Dny in range(1,2):#5): #range(Ny-10,Ny+1,2):
        k=Ny-Dny
        plt.clf()
        y1 = np.sin(k*x1)
        y2 = np.sin(k*x2)
        if 0:
            #plot the originals showing the beat frequency
            plt.plot(x1,y1,c='r')
            plt.plot(x2,y2,c='g')
        if 0:
            #plot the weird ass result
            other_nyquist=np.sin((Dny)*x1)
            super_nyquist=np.sin(Ny*x1)*np.cos(Dny*x1)-np.cos(Ny*x1)*np.sin(Dny*x1)
            plt.plot(x2,y2,c='g')
            plt.plot(x1,super_nyquist,c='k')
            #plt.plot(x2,y2-super_nyquist[::10],c='k')
        if 1:
            #plot the approximation
            #super_nyquist=np.sin(Ny*x1)*np.cos(Dny*x1)-np.cos(Ny*x1)*np.sin(Dny*x1)
            super_nyquist=-np.cos(Ny*x1)*np.sin(Dny*x1)
            plt.plot(x2,y2-super_nyquist[::10],c='k')
        #plt.plot(x1,other_nyquist,c='r')

        #plt.plot(x1,np.cos(Dny*x1),c='k')
        #plt.plot(x1,np.cos((Ny+Dny)*x1),c=[0.5]*3)
        #plt.plot(x1,,c=[0.5]*3)

        plt.savefig('%s/square_nyquist.pdf'%output_dir)
        time.sleep(1)

