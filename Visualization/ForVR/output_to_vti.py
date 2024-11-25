import numpy as np
from pyevtk.hl import imageToVTK

# Must install Pyevtk (https://github.com/paulo-herrera/PyEVTK) first:
#   `pip install pyevtk`

# References
#   VTK conversion: https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python
#   C_CONTIGUOUS flag: https://numpy.org/doc/stable/reference/generated/numpy.ascontiguousarray.html

def Output_scalar_field_to_vti( covering_grid, field_name, TakeLog=False, output_filename_prefix='', output_filename_suffix='', OrderZYX=False ):
    data_array        = np.array( covering_grid[field_name] )
    output_field_name = field_name[1] if isinstance(field_name, tuple) else field_name

    if TakeLog:
        data_array        = np.log10( data_array )
        output_field_name = 'Log_'+output_field_name

    if OrderZYX:
        # swap axes order from zyx to xyz in .vti file
        # make the flag C_CONTIGUOUS becomes True, to pass the vtk conversion assertion test
        data_array        = np.ascontiguousarray( np.swapaxes( data_array, 0, 2 ) )

    # NVIDIA Index only accepts point data for VTK files; cell data does not work for NVIDIA Index but works for Volume
    imageToVTK( '%s%s%s'%( output_filename_prefix, np.char.capitalize(output_field_name), output_filename_suffix ),
                pointData={ output_field_name : data_array.astype(np.float32) } )

def Output_vector_field_to_vti( covering_grid, vector_field, output_filename_prefix='', output_filename_suffix='', OrderZYX=False ):
    data_array_x      = np.array( covering_grid[vector_field+'_x'] )
    data_array_y      = np.array( covering_grid[vector_field+'_y'] )
    data_array_z      = np.array( covering_grid[vector_field+'_z'] )

    if OrderZYX:
        data_array_x  = np.ascontiguousarray( np.swapaxes( data_array_x, 0, 2 ) )
        data_array_y  = np.ascontiguousarray( np.swapaxes( data_array_y, 0, 2 ) )
        data_array_z  = np.ascontiguousarray( np.swapaxes( data_array_z, 0, 2 ) )

    imageToVTK( '%s%s%s'%( output_filename_prefix, np.char.capitalize(vector_field), output_filename_suffix ),
                pointData={ vector_field : ( data_array_x.astype(np.float32),
                                             data_array_y.astype(np.float32),
                                             data_array_z.astype(np.float32) ) } )

def Output_coordinate_field_to_vti( covering_grid, output_filename_prefix='', output_filename_suffix='', OrderZYX=False ):
    data_array_x      = np.array( covering_grid[('gamer','x')].in_units('code_length').d )
    data_array_y      = np.array( covering_grid[('gamer','y')].in_units('code_length').d )
    data_array_z      = np.array( covering_grid[('gamer','z')].in_units('code_length').d )

    if OrderZYX:
        data_array_x  = np.ascontiguousarray( np.swapaxes( data_array_x, 0, 2 ) )
        data_array_y  = np.ascontiguousarray( np.swapaxes( data_array_y, 0, 2 ) )
        data_array_z  = np.ascontiguousarray( np.swapaxes( data_array_z, 0, 2 ) )

    imageToVTK( '%sCoordX%s'%( output_filename_prefix, output_filename_suffix ),
                pointData={ 'Coord_x' : data_array_x.astype(np.float32) } )

    imageToVTK( '%sCoordY%s'%( output_filename_prefix, output_filename_suffix ),
                pointData={ 'Coord_y' : data_array_y.astype(np.float32) } )

    imageToVTK( '%sCoordZ%s'%( output_filename_prefix, output_filename_suffix ),
                pointData={ 'Coord_z' : data_array_z.astype(np.float32) } )

    imageToVTK( '%sVectorX%s'%( output_filename_prefix, output_filename_suffix ),
                pointData={ 'Vector_x' : ( 1.0*data_array_x.astype(np.float32),
                                           0.0*data_array_y.astype(np.float32),
                                           0.0*data_array_z.astype(np.float32) ) } )

    imageToVTK( '%sVectorY%s'%( output_filename_prefix, output_filename_suffix ),
                pointData={ 'Vector_y' : ( 0.0*data_array_x.astype(np.float32),
                                           1.0*data_array_y.astype(np.float32),
                                           0.0*data_array_z.astype(np.float32) ) } )

    imageToVTK( '%sVectorZ%s'%( output_filename_prefix, output_filename_suffix ),
                pointData={ 'Vector_z' : ( 0.0*data_array_x.astype(np.float32),
                                           0.0*data_array_y.astype(np.float32),
                                           1.0*data_array_z.astype(np.float32) ) } )
