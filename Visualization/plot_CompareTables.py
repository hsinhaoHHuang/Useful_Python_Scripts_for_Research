#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

####################################################
# Genereal Setting
####################################################
dpi        = 150
FONT_SIZE  = 16.0
LINE_WIDTH = 1.5
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
    file_table  = np.genfromtxt( filename, delimiter=None, comments='#', names=file_header, dtype=None, encoding=None )

    print( 'done!' )
    return file_table

def main() -> None:

    Table1 = LoadTable( '' )
    Table2 = LoadTable( '' )

    # create figure
    fig = plt.figure( figsize=(6,6) )

    # plot the lines with the colors
    ax = fig.add_subplot(111)

    ax.set_title( 'Record__MemInfo' )
    ax.plot( Table1['x'], Table1['y'], '.-', color='C1', label='Table1' )
    ax.plot( Table2['x'], Table2['y'], '.-', color='C2', label='Table2' )
    ax.set_xlabel( 'x' )
    ax.set_ylabel( 'y' )
    ax.legend()

    # save to file
    plt.tight_layout( pad=0.2 )
    fig.set_dpi(dpi)
    fig.savefig( 'fig_CompareTable.png', dpi=dpi )

if __name__ == '__main__':
    main()
