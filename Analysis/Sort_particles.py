
import numpy as np
import matplotlib.pyplot as plt

dtype_par = np.double

print( '' )
print( 'Loading %s ... '%'Record__Dump' )

with open( 'Record__Dump', 'r' ) as f:
    for line in f:
        if line.startswith( '#' ):
            header = line
        elif line.startswith( '\n' ):
            continue
        else:
            break #stop when there are no more #

    f.close()

RecordDump_header = header[1:].strip().split()
RecordDump_table  = np.genfromtxt( 'Record__Dump', delimiter=None, comments='#',
                                   names=RecordDump_header, dtype=None, encoding=None )

print( 'done!' )

for ID in RecordDump_table['DumpID']:

    Particle_Filename = 'Particle_%06d.txt'%(ID)
    print( '' )
    print('ID = %d'%(ID) )

    print( '' )
    print( 'Loading %s ... '%Particle_Filename )

    with open( Particle_Filename, 'r' ) as f:
        for line in f:
            if line.startswith( '#' ):
                header = line
            elif line.startswith( '\n' ):
                continue
            else:
                break #stop when there are no more #

        f.close()

    Particle_header = header[1:].strip().split()
    Particle_table  = np.genfromtxt( Particle_Filename, delimiter=None, comments='#',
                                     names=Particle_header, dtype=dtype_par, encoding=None )

    print( 'done!' )

    # sorting
    p1 = np.lexsort( (Particle_table['ParPosZ'], Particle_table['ParPosY'], Particle_table['ParPosX']) )

    # save the profile to text file
    np.savetxt( '%s_sorted'%(Particle_Filename),
                np.column_stack( (Particle_table['ParMass'][p1],
                                  Particle_table['ParPosX'][p1],
                                  Particle_table['ParPosY'][p1],
                                  Particle_table['ParPosZ'][p1],
                                  Particle_table['ParVelX'][p1],
                                  Particle_table['ParVelY'][p1],
                                  Particle_table['ParVelZ'][p1],
                                  Particle_table['ParTime'][p1] ) ),
                fmt='          % .16e',
                header='                        ParMass                           ParPosX                           ParPosY                           ParPosZ                           ParVelX                           ParVelY                           ParVelZ                           ParTime' )
