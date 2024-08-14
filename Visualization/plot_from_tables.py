#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

####################################################
# Genereal Setting
####################################################
dpi        = 150
FONT_SIZE  = 12.0
LINE_WIDTH = 2.0
MARKER_SIZE= 6.0

plt.rcParams['font.size']         = FONT_SIZE
plt.rcParams['figure.titlesize']  = FONT_SIZE
plt.rcParams['figure.dpi']        = dpi

plt.rcParams['axes.titlesize']    = FONT_SIZE
plt.rcParams['axes.labelsize']    = FONT_SIZE
plt.rcParams['axes.labelpad']     = FONT_SIZE/8.0*0.05
plt.rcParams['axes.linewidth']    = LINE_WIDTH

plt.rcParams['legend.fontsize']   = FONT_SIZE
plt.rcParams['lines.linewidth']   = LINE_WIDTH

plt.rcParams['xtick.major.size']  = 4.0*LINE_WIDTH
plt.rcParams['xtick.major.width'] = LINE_WIDTH
plt.rcParams['xtick.minor.size']  = 2.0*LINE_WIDTH
plt.rcParams['xtick.minor.width'] = 0.5*LINE_WIDTH
plt.rcParams['xtick.labelsize']   = FONT_SIZE

plt.rcParams['ytick.major.size']  = 4.0*LINE_WIDTH
plt.rcParams['ytick.major.width'] = LINE_WIDTH
plt.rcParams['ytick.minor.size']  = 2.0*LINE_WIDTH
plt.rcParams['ytick.minor.width'] = 0.5*LINE_WIDTH
plt.rcParams['ytick.labelsize']   = FONT_SIZE

plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'
###################################################################################################

def LoadTable( filename ):
    print( '' )
    print( 'Loading %s ... '%filename )

    with open( filename, 'r' ) as f:
        for line in f:
            if line.startswith( '#' ):
                header = line
            elif line.startswith( '\n' ):
                continue
            else:
                break #stop when there are no more #

        f.close()

    file_header = header[1:].strip().split()
    file_table  = np.genfromtxt( filename, delimiter=None, comments='#',
                                 names=file_header, dtype=None, encoding=None )

    print( 'done!' )
    return file_table

def main() -> None:

    # load data from table
    Table_1 = LoadTable( '/path/File_1' )
    Array_X = Table_1['X']

    Table_2 = LoadTable( '/path/File_2' )
    Array_Y = Table_2['Y']

    # create figure
    fig = plt.figure( figsize=(8,6) )

    ax = fig.add_subplot(111)

    ax.plot( Array_X, Array_Y, linestyle='-',  color='C0', label='XXX' )

    ax.set_title( 'XXX' )
    ax.set_xlabel( 'X' )
    ax.set_ylabel( 'Y' )
    ax.legend()

    #ax.set_xlim(,)
    #ax.set_ylim(,)
    ax.grid()

    ax.xaxis.get_ticklocs( minor=True )
    ax.yaxis.get_ticklocs( minor=True )
    ax.minorticks_on()
    ax.xaxis.set_ticks_position( 'both' )
    ax.yaxis.set_ticks_position( 'both' )
    ax.tick_params( which='both', direction='in' )

    # log scale
    #ax.set_xscale('log')
    #ax.set_yscale('log')

    # save to file
    plt.tight_layout( pad=0.2 )
    fig.set_dpi(dpi)
    fig.savefig( 'fig_from_tables.png', dpi=dpi )


if __name__ == '__main__':
    main()
