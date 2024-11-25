import argparse
import sys
import yt
import numpy as np
import add_ELBDM_derived_fields
import add_VR_derived_fields
import output_to_vti


# load the command-line parameters
parser = argparse.ArgumentParser( description='Extract uniform grid data from GAMER dataset\n\
                                               Example usage: python3 extract_uniform.py -s 0 -e 0 -p ../ -lv 0 -xL 0.0 -yL 0.0 -zL 0.0 -xR 1.0 -yR 1.0 -zR 1.0' )

parser.add_argument( '-s',  action='store', required=True,  type=int,   dest='idx_start', help='first data index'                                  )
parser.add_argument( '-e',  action='store', required=True,  type=int,   dest='idx_end',   help='last data index'                                   )
parser.add_argument( '-d',  action='store', required=False, type=int,   dest='didx',      help='delta data index [%(default)d]',      default=1    )
parser.add_argument( '-p',  action='store', required=False, type=str,   dest='prefix',    help='data path prefix [%(default)s]',      default='./' )
parser.add_argument( '-lv', action='store', required=False, type=int,   dest='lv',        help='sampling level [%(default)d]',        default=0    )
parser.add_argument( '-xL', action='store', required=False, type=float, dest='xL',        help='starting x coordinate [%(default)f]', default=-1.0 )
parser.add_argument( '-yL', action='store', required=False, type=float, dest='yL',        help='starting y coordinate [%(default)f]', default=-1.0 )
parser.add_argument( '-zL', action='store', required=False, type=float, dest='zL',        help='starting z coordinate [%(default)f]', default=-1.0 )
parser.add_argument( '-xR', action='store', required=False, type=float, dest='xR',        help='ending x coordinate [%(default)f]',   default=-1.0 )
parser.add_argument( '-yR', action='store', required=False, type=float, dest='yR',        help='ending y coordinate [%(default)f]',   default=-1.0 )
parser.add_argument( '-zR', action='store', required=False, type=float, dest='zR',        help='ending z coordinate [%(default)f]',   default=-1.0 )

args=parser.parse_args()


# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix
lv          = args.lv
xL          = args.xL
yL          = args.yL
zL          = args.zL
xR          = args.xR
yR          = args.yR
zR          = args.zR


# setting
PLOT_SLICES = True       # Whether to plot the slices
SMOOTHED    = True       # Whether to perform smoothing during extracting the covering grids
OUTPUT_H5   = True       # Whether to output hdf5 file
OUTPUT_NPZ  = True       # Whether to output npz file
OUTPUT_BIN  = True       # Whether to output binary file
OUTPUT_VTI  = True       # Whether to output .vti field (for VR)


# fields to be normalized by a radial profile
# - for example, to see the vortices in a halo
# - must be some kind of density fields
fields_prof_norm         = [ ('gamer', 'Dens'),
                             ('gamer', 'total_kinetic_energy_density'),
                             ('gamer', 'bulk_kinetic_energy_density'),
                             ('gamer', 'thermal_kinetic_energy_density'),
                           ]

# fields to be taken log before output as vti file
fields_output_log_to_vti = fields_prof_norm +\
                           [ ('gamer', 'profile_normalized_Dens'),
                             ('gamer', 'profile_normalized_total_kinetic_energy_density'),
                             ('gamer', 'profile_normalized_bulk_kinetic_energy_density'),
                             ('gamer', 'profile_normalized_thermal_kinetic_energy_density'),
                             ('gamer', 'reciprocal_profile_normalized_Dens'),
                           ]

# fields to extract
fields                   = [ ('gamer', 'total_velocity_magnitude'),
                             ('gamer', 'bulk_velocity_magnitude'),
                             ('gamer', 'thermal_velocity_magnitude'),
                             ('gamer', 'Real'),
                             ('gamer', 'Imag'),
                             ('gamer', 'S'),
                           ] + fields_output_log_to_vti

# fields to be saved as a vector filed in vti file
fields_output_vector_to_vti = [ 'bulk_velocity', 'thermal_velocity' ]


ts = yt.DatasetSeries( [ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

    # Get data information
    idx = int(('%s'%ds)[5:11])

    l_edge_x = ds.domain_left_edge[0]  if xL < 0.0 else ds.quan( xL, 'code_length' )
    l_edge_y = ds.domain_left_edge[1]  if yL < 0.0 else ds.quan( yL, 'code_length' )
    l_edge_z = ds.domain_left_edge[2]  if zL < 0.0 else ds.quan( zL, 'code_length' )
    r_edge_x = ds.domain_right_edge[0] if xR < 0.0 else ds.quan( xR, 'code_length' )
    r_edge_y = ds.domain_right_edge[1] if yR < 0.0 else ds.quan( yR, 'code_length' )
    r_edge_z = ds.domain_right_edge[2] if zR < 0.0 else ds.quan( zR, 'code_length' )

    N_x      = int(np.ceil( ds.domain_dimensions[0]*(2**lv) * (r_edge_x-l_edge_x)/ds.domain_width[0] ))
    N_y      = int(np.ceil( ds.domain_dimensions[1]*(2**lv) * (r_edge_y-l_edge_y)/ds.domain_width[1] ))
    N_z      = int(np.ceil( ds.domain_dimensions[2]*(2**lv) * (r_edge_z-l_edge_z)/ds.domain_width[2] ))

    r_edge_x = l_edge_x + N_x * ds.domain_width[0]/( ds.domain_dimensions[0]*(2**lv) )
    r_edge_y = l_edge_y + N_y * ds.domain_width[1]/( ds.domain_dimensions[1]*(2**lv) )
    r_edge_z = l_edge_z + N_z * ds.domain_width[2]/( ds.domain_dimensions[2]*(2**lv) )

    l_edge   = ds.arr( [ l_edge_x, l_edge_y, l_edge_z] )
    r_edge   = ds.arr( [ r_edge_x, r_edge_y, r_edge_z] )
    N        = ds.arr( [ N_x,       N_y,     N_z     ] )

    filename_prefix = 'CoveringGrid_%06d__'%idx
    filename_suffix = '__x%.3e--%.3e_y%.3e--%.3e_z%.3e--%.3e_Lv%02d_Nx%d_Ny%d_Nz%d'%( l_edge_x, r_edge_x, l_edge_y, r_edge_y, l_edge_z, r_edge_z, lv, N_x, N_y, N_z )


    # Add derived fields
    add_ELBDM_derived_fields.Add_ELBDM_derived_fields( ds )

    for field in fields_prof_norm:
        add_VR_derived_fields.Add_profile_normalized_field( ds, field, ds.domain_center, 0.4*ds.domain_width[0] )

    add_VR_derived_fields.Add_reciprocal_field( ds, ('gamer', 'profile_normalized_Dens') )


    # Plot slices
    if PLOT_SLICES:
        for field in fields:
            box = ds.box( l_edge, r_edge )
            box.max_level = lv
            sz = yt.SlicePlot( ds, 'z', field, center=0.5*( l_edge + r_edge ), data_source=box )
            sz.save( 'fig_Slice_%s.png'%field[1] )


    # Extract covering grid
    # reference: https://yt-project.org/docs/dev/examining/low_level_inspection.html#examining-grid-data-in-a-fixed-resolution-array
    if SMOOTHED:
        covering_grid = ds.smoothed_covering_grid( level=lv, left_edge=l_edge, dims=N, fields=fields, num_ghost_zones=16 )
    else:
        covering_grid =          ds.covering_grid( level=lv, left_edge=l_edge, dims=N, fields=fields, num_ghost_zones=16 )


    # Output data
    if OUTPUT_H5:
        covering_grid.save_as_dataset( filename='%sAll%s.h5'%( filename_prefix, filename_suffix ), fields=fields )
        # read later using h5py:
        #     file = h5py.File( 'filename.h5', 'r' )
        #     data = file['grid'][field]
        # or using yt:
        #     ds   = yt.load( 'filename.h5' )

    if OUTPUT_NPZ:
        for field in fields:
            # swap the axes order to save as an array with shape=(Nz, Ny, Nx) and the value can be accessed by arr[z][y][x]
            np.savez( '%s%s%s.npz'%( filename_prefix, np.char.capitalize(field[1]), filename_suffix ), Array=np.ascontiguousarray( np.swapaxes( covering_grid[field].astype(np.float32), 0, 2 ) ) )
            # read later using numpy: `arr = np.load( 'filename.npz' )['Array']`

    if OUTPUT_BIN:
        for field in fields:
            # swap the axes order to save as an array with shape=(Nz, Ny, Nx) and the value can be accessed by arr[z][y][x]
            np.ascontiguousarray( np.swapaxes( covering_grid[field].astype(np.float32), 0, 2 ) ).tofile( '%s%s%s.bin'%( filename_prefix, np.char.capitalize(field[1]), filename_suffix ) )
            # read later using numpy: `arr = np.fromfile( 'filename.bin', dtype=np.float32 ).reshape( Nz, Ny, Nx )`

    if OUTPUT_VTI:
        for field in fields:
            output_to_vti.Output_scalar_field_to_vti( covering_grid, field, TakeLog=False, output_filename_prefix=filename_prefix, output_filename_suffix=filename_suffix, OrderZYX=False )

        for field in fields_output_log_to_vti:
            output_to_vti.Output_scalar_field_to_vti( covering_grid, field, TakeLog=True, output_filename_prefix=filename_prefix, output_filename_suffix=filename_suffix, OrderZYX=False )

        for field in fields_output_vector_to_vti:
            output_to_vti.Output_vector_field_to_vti( covering_grid, field, filename_prefix, filename_suffix, OrderZYX=False )

        output_to_vti.Output_coordinate_field_to_vti( covering_grid, filename_prefix, filename_suffix, OrderZYX=False )
