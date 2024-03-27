#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np


def plot_example_figures_of_two_sets_of_time_evolution() -> None:

    # number of lines
    n_lines = 8

    # sampling points
    x    = np.linspace( 0.0, 4.0*np.pi,     100 )
    time = np.linspace( 0.0, 2.0*np.pi, n_lines )

    # create figure
    fig = plt.figure()
    ax  = fig.add_subplot(111)

    for i, t in enumerate(time):
       ax.plot( x, 0.1*t*np.sin(x  )+0.1*t, linestyle='-',  linewidth=2.0, color=plt.cm.autumn_r( 0.1+(i+1)*0.9/n_lines ), zorder=2.5, label='1, t = %.2f'%t )
       ax.plot( x, 0.1*t*np.sin(x-t)+0.1*t, linestyle='--', linewidth=2.5, color=plt.cm.winter_r( 0.1+(i+1)*0.9/n_lines ), zorder=2.0, label='2, t = %.2f'%t )

    ax.legend()

    ax.set_xlim( 0.0, 1.5*np.max(x) )

    # save to file
    plt.tight_layout()
    fig.savefig( 'fig_two_sets_of_time_evolution.png' )


def main() -> None:
    plot_example_figures_of_two_sets_of_time_evolution()

if __name__ == '__main__':
    main()
