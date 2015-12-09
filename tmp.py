
#norm = matplotlib.colors.Normalize(1,1e3)
plt.clf()
density_factor = 1000*4.6*3.08e18 #1000 cm^-3 * 4.6 pc * cm/pc
plt.imshow(factor1* np.log10( den), origin='lower', interpolation='nearest',cmap='gray')
cb=plt.colorbar()
cb.set_label(r'$N[\rm{cm}^{-2}]=N[\rm{code}]*%0.2e$'%density_factor)
dx = 4.6/8193.
old_xticks = plt.xticks()[0][1:-1]
plt.xticks( old_xticks, ["$%0.2f$"%n for n in  old_xticks*dx+L[0] ] )
plt.xlabel(r'$x[\rm{pc}]\  (\Delta x=$%s$\rm{pc})$'%expform(dx))
old_yticks = plt.yticks()[0][1:-1]
plt.yticks( old_yticks, ["$%0.2f$"%n for n in  old_yticks*dx+L[1] ] )
plt.ylabel(r'$y[\rm{pc}]\  (\Delta x=$%s$\rm{pc})$'%expform(dx))
#plt.yticks( plt.yticks()[0]*dx+L[1])


outname = 'p37_tick_test.png'
plt.savefig(outname)
print outname
