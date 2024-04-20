
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

    np.set_printoptions( formatter={'float': '{: .16e}'.format} )

    print( '---------------------------------------------------------------------' )
    print( '' )
    print( 'Particle' )
    print( 'NPar = ', Particle_table['ParMass'].size )
    print( 'Mass = ', Particle_table['ParMass'] )
    print( 'PosX = ', Particle_table['ParPosX'] )
    print( 'VelX = ', Particle_table['ParVelX'] )
    print( 'Time = ', Particle_table['ParTime'] )

    print( '' )
    print( 'Total Mass = {: .16e}'.format( np.sum( Particle_table['ParMass'] ) ) )
    print( 'Total MomX = {: .16e}'.format( np.sum( Particle_table['ParMass']*Particle_table['ParVelX'] ) ) )
    print( 'Total MomY = {: .16e}'.format( np.sum( Particle_table['ParMass']*Particle_table['ParVelY'] ) ) )
    print( 'Total MomZ = {: .16e}'.format( np.sum( Particle_table['ParMass']*Particle_table['ParVelZ'] ) ) )
