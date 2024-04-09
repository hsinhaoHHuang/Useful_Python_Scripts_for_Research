#!/usr/bin/env python3

from mpi4py import MPI
import numpy as np


################################################################################################################
# Send and Recv
################################################################################################################
# in real code, this section might
# read in data parameters from a file
'''
################################################################################################################
# Bcast
################################################################################################################
if rank == 0:
    # create a data array on process 0
    # in real code, this section might
    # read in data parameters from a file
    numData_bcasting = 10
    data_bcasting = np.linspace(0.0,3.14,numData_bcasting)
else:
    numData_bcasting = None

# broadcast numData and allocate array on other ranks:
numData_bcasting = comm.bcast(numData_bcasting, root=0)
if rank != 0:
    data_bcasting = np.empty(numData_bcasting, dtype='d')

comm.Bcast(data_bcasting, root=0) # broadcast the array from rank 0 to all others

print('Rank: ',rank, ', data received: ',data_bcasting)

################################################################################################################
# Scattering
################################################################################################################
numDataPerRank_scattering = 10
data_scattering = None
if rank == 0:
    data_scattering = np.linspace(1,size*numDataPerRank_scattering,numDataPerRank_scattering*size)
    # when size=4 (using -n 4), data = [1.0:40.0]

recvbuf = np.empty(numDataPerRank_scattering, dtype='d') # allocate space for recvbuf
comm.Scatter(data_scattering, recvbuf, root=0)

print('Rank: ',rank, ', recvbuf received: ',recvbuf)

################################################################################################################
# Gathering
################################################################################################################
numDataPerRank_gathering = 10
sendbuf = np.linspace(rank*numDataPerRank_gathering+1,(rank+1)*numDataPerRank_gathering,numDataPerRank_gathering)
print('Rank: ',rank, ', sendbuf: ',sendbuf)

recvbuf = None
if rank == 0:
    recvbuf = np.empty(numDataPerRank_gathering*size, dtype='d')

comm.Gather(sendbuf, recvbuf, root=0)

if rank == 0:
    print('Rank: ',rank, ', recvbuf received: ',recvbuf)

################################################################################################################
# Reduce
################################################################################################################
# Create some np arrays on each process:
# For this demo, the arrays have only one
# entry that is assigned to be the rank of the processor
value = np.array(rank,'d')

print(' Rank: ',rank, ' value = ', value)

# initialize the np arrays that will store the results:
value_sum      = np.array(0.0,'d')
value_max      = np.array(0.0,'d')

# perform the reductions:
comm.Reduce(value, value_sum, op=MPI.SUM, root=0)
comm.Reduce(value, value_max, op=MPI.MAX, root=0)

if rank == 0:
    print(' Rank 0: value_sum =    ',value_sum)
    print(' Rank 0: value_max =    ',value_max)
################################################################################################################
'''

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
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    print( f'My rank is Rank-{rank:02d}' )

    numData_sending = 10
    comm.send(numData_sending, dest=(rank+1)%size)

    data_sending = rank*np.linspace(0.0,3.14,numData_sending)
    comm.Send(data_sending, dest=(rank+1)%size)


numData_recving = comm.recv(source=(rank+size-1)%size)
print('Number of data to receive: ',numData_recving)

data_recving = np.empty(numData_recving, dtype='d')  # allocate space to receive the array
comm.Recv(data_recving, source=(rank+size-1)%size)

print('data received: ',data_recving)


    print('Hello World!')

if __name__ == '__main__':
    main()
