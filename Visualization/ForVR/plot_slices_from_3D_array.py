import matplotlib.pyplot as plt
import numpy as np

from matplotlib.colors import Normalize, LogNorm


def imshow3d(ax, array, value_direction='z', pos=0, norm=None, cmap=None):
    """
    Display a 2D array as a  color-coded 2D image embedded in 3d.

    The image will be in a plane perpendicular to the coordinate axis *value_direction*.

    Parameters
    ----------
    ax : Axes3D
        The 3D Axes to plot into.
    array : 2D numpy array
        The image values.
    value_direction : {'x', 'y', 'z'}
        The axis normal to the image plane.
    pos : float
        The numeric value on the *value_direction* axis at which the image plane is
        located.
    norm : `~matplotlib.colors.Normalize`, default: Normalize
        The normalization method used to scale scalar data. See `imshow()`.
    cmap : str or `~matplotlib.colors.Colormap`, default: :rc:`image.cmap`
        The Colormap instance or registered colormap name used to map scalar data
        to colors.
    """
    if norm is None:
        norm = Normalize()
    colors = plt.get_cmap(cmap)(norm(array))

    if value_direction == 'x':
        ny, nz = array.shape
        yi, zi = np.mgrid[0:ny + 1, 0:nz + 1]
        xi = np.full_like(yi, pos)
    elif value_direction == 'y':
        nz, nx = array.shape
        zi, xi = np.mgrid[0:nz + 1, 0:nx + 1]
        yi = np.full_like(zi, pos)
    elif value_direction == 'z':
        nx, ny = array.shape
        xi, yi = np.mgrid[0:nx + 1, 0:ny + 1]
        zi = np.full_like(xi, pos)
    else:
        raise ValueError(f"Invalid value_direction: {value_direction!r}")
    ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, facecolors=colors, shade=False)


# from .npz
filename = 'CoveringGrid.npz'
data = np.load(filename)['Array']
Nz, Ny, Nx = data.shape

# from .bin
#filename = 'CoveringGrid.bin'
#Nx =
#Ny =
#Nz =
#data = np.fromfile( filename, dtype=np.float32 ).reshape( Nz, Ny, Nx )

# from .h5
#filename  = 'CoveringGrid.h5'
#file_h5py = h5py.File( filename, 'r' )
#data = np.ascontiguousarray( np.swapaxes( np.array( file_h5py['grid']['Dens'].in_units('code_density').d ), 0, 2 ) )
#Nz, Ny, Nx = data.shape

X0 = Nx//2
Y0 = Ny//2
Z0 = Nz//2

data_xy = np.swapaxes( data[Z0, :,  :], 0, 1 )
data_yz = np.swapaxes( data[:,  :, X0], 0, 1 )
data_zx =              data[:,  Y0, :]

print( f'{data.shape    = }' )
print( f'{data_xy.shape = }' )
print( f'{data_yz.shape = }' )
print( f'{data_xz.shape = }' )

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set(xlabel="x", ylabel="y", zlabel="z")

imshow3d( ax, data_xy, value_direction='z', pos=Z0, cmap='viridis', norm=LogNorm(vmin=data.min(), vmax=data.max(), clip=True) )
imshow3d( ax, data_yz, value_direction='x', pos=X0, cmap='viridis', norm=LogNorm(vmin=data.min(), vmax=data.max(), clip=True) )
imshow3d( ax, data_zx, value_direction='y', pos=Y0, cmap='viridis', norm=LogNorm(vmin=data.min(), vmax=data.max(), clip=True) )

# Make the grid
#x, y, z = np.meshgrid(np.arange(-0.8, 1, 0.2),
#                      np.arange(-0.8, 1, 0.2),
#                      np.arange(-0.8, 1, 0.8))
#
## Make the direction data for the arrows
#u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
#v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
#w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
#     np.sin(np.pi * z))
#
#ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)



plt.show()
