
ef('p37_tools.py')
distance = 150. #pc
resolution = 40.#arcsec
smooth_pc = distance*resolution/206264. #in pc
smooth_px = max(int(smooth_pc/(4.6/Nzones)), 2)
ef('p37_tools.py')
blur = blur_image(den,smooth_px)
plt.clf()
if 1:
    plt.imshow( np.log10( den), origin='lower', interpolation='nearest')
    plt.savefig('full_a.png')
    plt.clf()
    plt.imshow( np.log10( blur), origin='lower', interpolation='nearest')
    plt.savefig('full_b.png')

