#!/usr/bin/env python3

from mpi4py import MPI
import numpy as np


def main() -> None:
    """main function
    This function is a template.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    # Initialization
    comm      = MPI.COMM_WORLD
    Nrank     = comm.Get_size()

    my_rank   = comm.Get_rank()
    prev_rank = (my_rank+Nrank-1)%Nrank
    next_rank = (my_rank      +1)%Nrank

    print( f'My rank is Rank-{my_rank:02d}' )

    # Send and Recv
    print( 'Send and Recv' )
    numData_sending = 10
    comm.send( numData_sending, dest=next_rank )

    data_sending = my_rank*np.linspace( 0.0, 3.14, numData_sending )
    comm.Send( data_sending, dest=next_rank )
    print( f'Rank-{my_rank:02d} data sent: ', data_sending )

    numData_recving = comm.recv( source=prev_rank )
    print( 'Number of data to receive: ', numData_recving )

    data_recving = np.empty( numData_recving, dtype='d' )  # allocate space to receive the array
    comm.Recv( data_recving, source=prev_rank )

    print( f'Rank-{my_rank:02d} data received: ', data_recving )

    # Bcast
    print( 'Bcase' )
    if my_rank == 0:
        # create a data array on process 0
        # in real code, this section might
        # read in data parameters from a file
        numData_bcasting = 10
        data_bcasting    = np.linspace( 0.0, 3.14, numData_bcasting )
    else:
        numData_bcasting = None

    # broadcast numData and allocate array on other my_ranks:
    numData_bcasting = comm.bcast( numData_bcasting, root=0 )

    if my_rank != 0:
        data_bcasting = np.empty( numData_bcasting, dtype='d' )

    comm.Bcast( data_bcasting, root=0 ) # broadcast the array from rank 0 to all others

    print( f'Rank-{my_rank:02d} data bcasted: ', data_bcasting )

    # Scattering
    print( 'Scatter' )
    numDataPerRank_scattering = 10
    data_scattering = None

    if my_rank == 0:
        data_scattering = np.linspace( 1, Nrank*numDataPerRank_scattering, Nrank*numDataPerRank_scattering )
        # when Nrank=4 (using -n 4), data = [1.0:40.0]

    recvbuf = np.empty( numDataPerRank_scattering, dtype='d' ) # allocate space for recvbuf
    comm.Scatter( data_scattering, recvbuf, root=0 )

    print( f'Rank-{my_rank:02d} data scattered: ', recvbuf )

    # Gathering
    print( 'Gather' )
    numDataPerRank_gathering = 10
    sendbuf = np.linspace( my_rank*numDataPerRank_gathering+1, (my_rank+1)*numDataPerRank_gathering, numDataPerRank_gathering )
    print( f'Rank-{my_rank:02d} data gathering: ', sendbuf )

    recvbuf = None
    if my_rank == 0:
        recvbuf = np.empty( numDataPerRank_gathering*Nrank, dtype='d' )

    comm.Gather( sendbuf, recvbuf, root=0 )

    if my_rank == 0:
        print( f'Rank-{my_rank:02d} data gathered: ', recvbuf )

    # Reduce
    print( 'Reduce' )
    value = np.array( my_rank, 'd' )

    print( ' Rank: ', my_rank, ' value = ', value )

    # initialize the np arrays that will store the results:
    value_sum   = np.array( 0.0, 'd' )
    value_max   = np.array( 0.0, 'd' )

    # perform the reductions:
    comm.Reduce( value, value_sum, op=MPI.SUM, root=0 )
    comm.Reduce( value, value_max, op=MPI.MAX, root=0 )

    if my_rank == 0:
        print( ' Rank 0: value_sum =    ', value_sum )
        print( ' Rank 0: value_max =    ', value_max )


if __name__ == '__main__':
    main()
