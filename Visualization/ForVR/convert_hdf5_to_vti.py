import argparse
import h5py
import numpy as np
import os
import output_to_vti

parser = argparse.ArgumentParser( description='Convert HDF5 into VTI \n\
                                  Example usage: python3 convert_hdf5_to_vti.py -input input.h5 -field Dens -log 1')

parser.add_argument( '-input', action='store', required=True,  type=str, dest='input_file',   help='input file name'                   )
parser.add_argument( '-field', action='store', required=True,  type=str, dest='target_field', help='target field name'                 )
parser.add_argument( '-log',   action='store', required=False, type=int, dest='take_log',     help='take log [%(default)d]', default=0 )

args=parser.parse_args()

input_file   = args.input_file
target_field = args.target_field
take_log     = args.take_log

file_h5py    = h5py.File( input_file, 'r' )
file_name    = os.path.splitext(os.path.basename(input_file))[0]

# from yt covering_grid.save_as_dataset()
grid_data    = file_h5py['grid']
order_zyx    = False

# from GAMER_ExtractUniform
# grid_data  = file_h5py['Data']
# order_zyx  = True

output_to_vti.Output_scalar_field_to_vti( grid_data, target_field, TakeLog=take_log, output_filename_prefix=file_name+'_', output_filename_suffix='', OrderZYX=order_zyx )
