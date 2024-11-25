import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d


filename      = 'CoveringGrid_000001__All__x3.235e-02--6.469e-02_y3.235e-02--6.469e-02_z3.235e-02--6.469e-02_Lv01_Nx320_Ny320_Nz320.h5'
#field         = ('grid','reciprocal_profile_normalized_Dens')
field         = ('grid','S')

ds = yt.load( filename )

dh = (ds.domain_width[0]/ds.domain_dimensions[0])
N = 1024

frb = yt.OffAxisSlicePlot( ds,
                           [ 0.670918, -0.232255, 0.70422 ],
                           field,
                           center=[ 320+102.206, 320+23.3333, 320+27.9842 ]*dh,
                           width=10*dh,
                           buff_size=[N, N] ).frb

fig = plt.figure(figsize=plt.figaspect(2.))
ax = fig.add_subplot(2,1,1,projection='3d')

X, Y = np.meshgrid( np.linspace( 0, 1, N ), np.linspace( 0, 1, N ) )

#Z = np.log10( frb[field].d )
Z = frb[field].d
z_max =  np.pi
z_min = -np.pi

# Plot the 3D surface
ax.plot_surface( X, Y, Z, cmap=cm.rainbow, vmin=z_min, vmax=z_max, linewidth=0, antialiased=False )

# Plot projections of the contours for each dimension.  By choosing offsets
# that match the appropriate axes limits, the projected contours will sit on
# the 'walls' of the graph
#ax.contourf(X, Y, Z, zdir='z', offset=2*z_max, cmap='rainbow')
#ax.contourf(X, Y, Z, zdir='x', offset=0, cmap='coolwarm')
#ax.contourf(X, Y, Z, zdir='y', offset=1, cmap='coolwarm')
#ax.quiver(X, Y, np.cos(X), np.sin(Y), units='width')

#ax.set( xlim=(0, 1), ylim=(0, 1), zlim=(z_min, z_max), xlabel='X', ylabel='Y', zlabel=r'log(1/$\rho)$' )
ax.set( xlim=(0, 1), ylim=(0, 1), zlim=(z_min, z_max), xlabel='X', ylabel='Y', zlabel='S' )

ax2 = fig.add_subplot(2,1,2)
im = ax2.imshow( Z, cmap='rainbow', vmin=z_min, vmax=z_max, interpolation='none', origin='lower' )
fig.colorbar(im, ax=ax2)
#ax2.quiver( x, y, np.cos(x), np.sin(y), scale=1, width=0.015 )

plt.show()
