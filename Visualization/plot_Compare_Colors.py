#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

Miville_color = {   'blue':'#1E90FF',
                     'red':'#FF6346',
                   'green':'#00FF01',
                   'black':'#000000',
                  'yellow':'#FFD600' }
Sasha_color = {      'red':'#e6194B',
                   'green':'#3cb44b',
                  'yellow':'#ffe119',
                    'blue':'#4363d8',
                  'orange':'#f58231',
                  'purple':'#911eb4',
                    'cyan':'#42d4f4',
                 'magenta':'#f032e6',
                    'lime':'#bfef45',
                    'pink':'#fabed4',
                    'teal':'#469990',
                 'lavende':'#dcbeff',
                   'brown':'#9A6324',
                   'beige':'#fffac8',
                  'maroon':'#800000',
                    'mint':'#aaffc3',
                   'olive':'#808000',
                 'apricot':'#ffd8b1',
                    'navy':'#000075',
                    'grey':'#a9a9a9',
                   'white':'#ffffff',
                   'black':'#000000' }

def plot_compare_colors() -> None:

    # sampling points
    x = np.linspace(0.0, 2.0*np.pi, 100)

    # create figure
    fig = plt.figure(figsize=(15,12))

    # plot the lines with the colors
    ax1 = fig.add_subplot(541)
    ax1.plot( x, np.sin(x+0), color='r',                  label='default' )
    ax1.plot( x, np.sin(x+1), color='xkcd:red',           label='xkcd'    )
    ax1.plot( x, np.sin(x+2), color=Miville_color['red'], label='Miville' )
    ax1.plot( x, np.sin(x+3), color=Sasha_color['red'],   label='Sasha'   )
    ax1.legend()

    ax2 = fig.add_subplot(542)
    ax2.plot( x, np.sin(x+0), color='g',                    label='default' )
    ax2.plot( x, np.sin(x+1), color='xkcd:green',           label='xkcd'    )
    ax2.plot( x, np.sin(x+2), color=Miville_color['green'], label='Miville' )
    ax2.plot( x, np.sin(x+3), color=Sasha_color['green'],   label='Sasha'   )
    ax2.legend()

    ax3 = fig.add_subplot(543)
    ax3.plot( x, np.sin(x+0), color='b',                    label='default' )
    ax3.plot( x, np.sin(x+1), color='xkcd:blue',            label='xkcd'    )
    ax3.plot( x, np.sin(x+2), color=Miville_color['blue'],  label='Miville' )
    ax3.plot( x, np.sin(x+3), color=Sasha_color['blue'],    label='Sasha'   )
    ax3.legend()

    ax4 = fig.add_subplot(544)
    ax4.plot( x, np.sin(x+0), color='y',                    label='default' )
    ax4.plot( x, np.sin(x+1), color='xkcd:yellow',          label='xkcd'    )
    ax4.plot( x, np.sin(x+2), color=Miville_color['yellow'],label='Miville' )
    ax4.plot( x, np.sin(x+3), color=Sasha_color['yellow'],  label='Sasha'   )
    ax4.legend()

    # save to file
    plt.tight_layout()
    fig.savefig( 'fig_compare_colors.png' )


def main() -> None:
    plot_compare_colors()

if __name__ == '__main__':
    main()
