#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np


def plot_example_figures_with_gradient_color_lines() -> None:

    # sampling points
    x = np.linspace(0.0, 2.0*np.pi, 100)

    # number of lines
    n_lines = 20

    # colors from different colormap
    colors_rainbow = plt.cm.rainbow( np.linspace(0, 1, n_lines) )
    colors_jet     = plt.cm.jet(     np.linspace(0, 1, n_lines) )
    colors_viridis = plt.cm.viridis( np.linspace(0, 1, n_lines) )
    colors_inferno = plt.cm.inferno( np.linspace(0, 1, n_lines) )
    colors_hsv     = plt.cm.hsv(     np.linspace(0, 1, n_lines) )
    colors_RdBu    = plt.cm.RdBu(    np.linspace(0, 1, n_lines) )

    # create figure
    fig = plt.figure(figsize=(12,6))

    # plot the lines with the colors
    ax1 = fig.add_subplot(321)
    ax1.set_title('rainbow')
    for i, c in enumerate(colors_rainbow):
       ax1.plot(x, np.sin(x)+0.1*i, color=c)

    ax2 = fig.add_subplot(322)
    ax2.set_title('jet')
    for i, c in enumerate(colors_jet):
       ax2.plot(x, np.sin(x)+0.1*i, color=c)

    ax3 = fig.add_subplot(323)
    ax3.set_title('viridis')
    for i, c in enumerate(colors_viridis):
       ax3.plot(x, np.sin(x)+0.1*i, color=c)

    ax4 = fig.add_subplot(324)
    ax4.set_title('inferno')
    for i, c in enumerate(colors_inferno):
       ax4.plot(x, np.sin(x)+0.1*i, color=c)

    ax5 = fig.add_subplot(325)
    ax5.set_title('hsv')
    for i, c in enumerate(colors_hsv):
       ax5.plot(x, np.sin(x)+0.1*i, color=c)

    ax6 = fig.add_subplot(326)
    ax6.set_title('RdBu')
    for i, c in enumerate(colors_RdBu):
       ax6.plot(x, np.sin(x)+0.1*i, color=c)

    # save to file
    plt.tight_layout()
    fig.savefig( 'fig_gradient_color_lines.png' )


def main() -> None:
    plot_example_figures_with_gradient_color_lines()

if __name__ == '__main__':
    main()
