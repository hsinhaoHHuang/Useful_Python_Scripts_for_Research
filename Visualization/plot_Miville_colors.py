#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['axes.prop_cycle'] = plt.cycler(color=[ '#1E90FF',
                                                     '#FF6346',
                                                     '#00FF01',
                                                     '#000000',
                                                     '#FFD600',
                                                   ])

def plot_example_Miville_colors() -> None:

    # sampling points
    x = np.linspace(0.0, 2.0*np.pi, 100)

    # number of lines
    n_lines = 10

    # create figure
    fig = plt.figure(figsize=(4,4))

    # plot the lines with the colors
    ax1 = fig.add_subplot(111)
    for i in range(10):
       ax1.plot( x, np.sin(x)+0.1*i )

    # save to file
    plt.tight_layout()
    fig.savefig( 'fig_Miville_colors.png' )


def main() -> None:
    plot_example_Miville_colors()

if __name__ == '__main__':
    main()
