#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np



def main() -> None:
    """main function
    This function plots a example meshgrid.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    # sampling points
    x = np.linspace(0.0, 1.0, 100)
    y = np.linspace(0.0, 1.0, 100)

    xx, yy = np.meshgrid(x, y)
    zz = np.sqrt(xx**2 + yy**2)

    # create figure
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    pos = ax1.imshow(zz, origin='lower', extent=[0,1,0,1], interpolation='none' )

    ax2 = fig.add_subplot(122)
    ax2.contourf(zz)
    ax2.axis('scaled')

    # save to file
    plt.tight_layout()
    fig.savefig( 'fig_meshgrid.png' )

if __name__ == '__main__':
    main()
