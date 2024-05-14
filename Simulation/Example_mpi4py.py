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

    comm.Barrier()
    # Send and Recv
    if my_rank == 0:
        print( '' )
        print( 'Send and Recv' )

    numData_sending = 10
    comm.send( numData_sending, dest=next_rank )

    data_sending    = my_rank*np.ones( numData_sending )
    comm.Send( data_sending, dest=next_rank )
    print( f'Rank-{my_rank:02d} data sent:     ', data_sending )

    numData_recving = comm.recv( source=prev_rank )
    data_recving    = np.empty( numData_recving, dtype='d' )  # allocate space to receive the array
    comm.Recv( data_recving, source=prev_rank )
    print( f'Rank-{my_rank:02d} data received: ', data_recving )

    comm.Barrier()
    # Bcast
    if my_rank == 0:
        print( '' )
        print( 'Bcast' )
    if my_rank == 0:
        # create a data array on process 0
        # in real code, this section might
        # read in data parameters from a file
        numData_bcasting = 10
        data_bcasting    = np.linspace( 0.0, 0.9, numData_bcasting )
        print( f'Rank-{my_rank:02d} data bcasting: ', data_bcasting )
    else:
        numData_bcasting = None

    # broadcast numData and allocate array on other my_ranks:
    numData_bcasted = comm.bcast( numData_bcasting, root=0 )

    if my_rank != 0:
        data_bcasting = np.empty( numData_bcasted, dtype='d' )

    comm.Bcast( data_bcasting, root=0 ) # broadcast the array from rank 0 to all others

    print( f'Rank-{my_rank:02d} data bcasted:  ', data_bcasting )

    comm.Barrier()
    # Scattering
    if my_rank == 0:
        print( '' )
        print( 'Scatter' )

    numDataPerRank_scattering = 10

    if my_rank == 0:
        data_scattering = np.linspace( 1, Nrank*numDataPerRank_scattering, Nrank*numDataPerRank_scattering )
        print( f'Rank-{my_rank:02d} data scattering: ', data_scattering )
        # when Nrank=4 (using -n 4), data = [1.0:40.0]
    else:
        data_scattering = None

    data_scattered = np.empty( numDataPerRank_scattering, dtype='d' ) # allocate space

    comm.Scatter( data_scattering, data_scattered, root=0 )
    print( f'Rank-{my_rank:02d} data scattered:  ', data_scattered )

    comm.Barrier()
    # Gathering
    if my_rank == 0:
        print( '' )
        print( 'Gather' )
    numDataPerRank_gathering = 10
    data_gathering = np.linspace( my_rank*numDataPerRank_gathering+1, (my_rank+1)*numDataPerRank_gathering, numDataPerRank_gathering )
    print( f'Rank-{my_rank:02d} data gathering: ', data_gathering )

    data_gathered = None
    if my_rank == 0:
        data_gathered = np.empty( numDataPerRank_gathering*Nrank, dtype='d' )

    comm.Gather( data_gathering, data_gathered, root=0 )

    if my_rank == 0:
        print( f'Rank-{my_rank:02d} data gathered:  ', data_gathered )

    comm.Barrier()
    # Reduce
    if my_rank == 0:
        print( '' )
        print( 'Reduce' )

    data_reducing = np.array( my_rank, 'd' )

    print( f'Rank-{my_rank:02d} data_reducing = ', data_reducing )

    # initialize the np arrays that will store the results:
    value_sum   = np.array( 0.0, 'd' )
    value_max   = np.array( 0.0, 'd' )

    # perform the reductions:
    comm.Reduce( data_reducing, value_sum, op=MPI.SUM, root=0 )
    comm.Reduce( data_reducing, value_max, op=MPI.MAX, root=0 )

    if my_rank == 0:
        print( ' Rank 0: value_sum =    ', value_sum )
        print( ' Rank 0: value_max =    ', value_max )


if __name__ == '__main__':
    main()
