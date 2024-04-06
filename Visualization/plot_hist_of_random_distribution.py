#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    """
    This function plots random points following the 2D normal distribution.

    Parameters:
        None

    Return:
        None
    """
    # create the figure
    fig = plt.figure()
    ax  = fig.add_subplot(111)

    lambda_poisson = 5
    num_points     = 10000

    Data = np.random.poisson( lambda_poisson, num_points )

    # plot the histogram
    ax.hist( Data, bins=range(3*lambda_poisson), alpha=0.7, color='blue', align='left' )

    # annotate the information
    text = ax.text( 2*lambda_poisson, 1700, 'N = %d\nlambda = %4.1f'%(num_points, lambda_poisson), color='k' )

    # set lables
    ax.set_xlabel( 'Number of events' )
    ax.set_ylabel( 'Frequency' )

    # save the figure
    plt.tight_layout()
    fig.savefig( 'fig_hist_of_random_distribution.png' )

if __name__ == '__main__':
    main()
