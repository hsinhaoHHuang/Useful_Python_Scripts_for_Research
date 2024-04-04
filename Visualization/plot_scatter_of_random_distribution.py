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

    mean       = 50.0
    sigma      = 10.0
    num_points = 10000

    X = np.random.normal( mean, sigma, num_points )
    Y = np.random.normal( mean, sigma, num_points )

    # plot the scatter
    ax.scatter( X, Y, 50, '0.0', lw=2            )
    ax.scatter( X, Y, 50, '1.0', lw=0            )
    ax.scatter( X, Y, 40, 'C1',  lw=0, alpha=0.1 )

    # annotate the information
    text = ax.text( mean+0.5*sigma, mean-0.5*sigma, 'mean = %4.3f\nsigma = %4.3f'%(mean, sigma), color='C0' )

    # set lables
    ax.set_xlabel( r'$x$' )
    ax.set_ylabel( r'$y$' )

    # x,y limit
    ax.set_xlim( mean-3*sigma, mean+3*sigma )
    ax.set_ylim( mean-3*sigma, mean+3*sigma )

    # save the figure
    plt.tight_layout()
    fig.savefig( 'fig_scatter_of_random_distribution.png' )

if __name__ == '__main__':
    main()
