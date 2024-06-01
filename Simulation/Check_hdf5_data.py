import h5py
import numpy as np


# input filename
f = h5py.File( 'Data_000000', 'r' )


for group in list(f):
    print( '- Group ', group )

    for dset in list(f[group]):
        print( '           - dset ', dset )

        for attribute in list( f[group][dset].attrs ):
            print( '                     - attribute ', attribute, '=', f[group][dset].attrs[attribute] )

        if group == 'GridData':
            GID = 0
            k   = 0
            j   = 0
            i   = 0
            print( '                     - size   = ', len( list( f[group][dset] ) ) )
            print( '                     - value  = ', f[group][dset][GID][k][j][i], ' at GID =', GID, 'k =', k, 'j =', j, 'i =', i )

        elif group == 'Particle':
            print( '                     - size   = ', len( list( f[group][dset] ) ) )
            print( '                     - value  = ', f[group][dset][0:10].tolist() )

        elif group == 'Info':
            for name in list(f[group][dset].dtype.names):
               print( '                     - %-30s ='%name, np.array( f[group][dset][name] ).tolist() )

        elif group == 'Tree':
            print( '                     - shape    =', f[group][dset].shape       )
            print( '                     - value[0] =', f[group][dset][0].tolist() )

        print( '' )

    print( '' )
