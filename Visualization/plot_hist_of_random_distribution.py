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
    ax.hist( Data, bins=range(15), alpha=0.7, color='blue', align='left' )

    # annotate the information
    #text = ax.text( mean+0.5*sigma, mean-0.5*sigma, 'mean = %4.3f\nsigma = %4.3f'%(mean, sigma), color='C0' )

    # set lables
    ax.set_xlabel( r'$Number of events$' )
    ax.set_ylabel( r'$Frequency$' )

    # x,y limit
    #ax.set_xlim( mean-3*sigma, mean+3*sigma )
    #ax.set_ylim( mean-3*sigma, mean+3*sigma )

    # save the figure
    plt.tight_layout()
    fig.savefig( 'fig_hist_of_random_distribution.png' )

if __name__ == '__main__':
    main()
