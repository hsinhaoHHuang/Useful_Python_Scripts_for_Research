#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, NullFormatter



####################################################
# Genereal Setting
####################################################
dpi        = 600
FONT_SIZE  = 8.0
LINE_WIDTH = 1.0
MARKER_SIZE= 4.0

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


####################################################
# Setting for Plotting
####################################################
n_rows        = 2           # number of rows in the figure panels
n_cols        = 3           # number of columns in the figure panels

color_1       = 'C3'        # color for Table_1
color_2       = 'C0'        # color for Table_2
color_ref     = 'C2'        # color for Reference

axis_y_min    = 8.0e+2      # min for y axis
axis_y_max    = 2.0e+6      # max for y axis
axis_x_min    = 5.0e-5      # min for x axis
axis_x_max    = 2.0e-2      # max for x axis

axis_x_label  = r'$\mathit{r}\ ({\rm kpc})}$'
axis_y_label  = r'$\mathit{\rho}\ ({\rm M}_{\odot}{\rm kpc}^{\rm -3})$'

annotated_arrow_1_x = 1.0
annotated_arrow_1_y = 1.0
annotated_text_1_x  = 1.0
annotated_text_1_y  = 1.0
annotated_text_1    = '1'

annotated_arrow_2_x = 2.0
annotated_arrow_2_y = 2.0
annotated_text_2_x  = 2.0
annotated_text_2_y  = 2.0
annotated_text_2    = '2'

cm            = 1/2.54      # centimeters in inches
fig_size_x    = 16*cm       # output figure size
fig_size_y    = (fig_size_x-1*cm)/n_cols*n_rows
fig           = plt.figure( 1, (fig_size_x, fig_size_y), dpi=dpi )
###################################################################################################


####################################################
# Read the data and parameters
####################################################
# path
Table_1_Path     = np.empty( n_rows*n_cols, dtype=object )
Table_1_Path[0]  = './Table_1_000000'
Table_1_Path[1]  = './Table_1_000001'
Table_1_Path[2]  = './Table_1_000002'
Table_1_Path[3]  = './Table_1_000003'
Table_1_Path[4]  = './Table_1_000004'
Table_1_Path[5]  = './Table_1_000005'

Table_2_Path     = np.empty( n_rows*n_cols, dtype=object )
Table_2_Path[0]  = './Table_2_000000'
Table_2_Path[1]  = './Table_2_000001'
Table_2_Path[2]  = './Table_2_000002'
Table_2_Path[3]  = './Table_2_000003'
Table_2_Path[4]  = './Table_2_000004'
Table_2_Path[5]  = './Table_2_000005'
###################################################################################################


#  Define Functions
def Analytical_Ref(x):

    y_analytical = 1e1*x**(-1)

    return y_analytical



def main() -> None:

    # Loop for each panel
    for panel_idx in range(0, n_rows*n_cols, 1):

        # load the table
        Table_1_Data_x, Table_1_Data_y = np.loadtxt( Table_1_Path[panel_idx], unpack=True )
        Table_2_Data_x, Table_2_Data_y = np.loadtxt( Table_2_Path[panel_idx], unpack=True )

        ax  = fig.add_subplot( n_rows, n_cols, panel_idx+1 )

        # plot data
        ax.plot( Table_1_Data_x, Table_1_Data_y,             linestyle='-',  color=color_1,    linewidth=LINE_WIDTH, label='1' )
        ax.plot( Table_2_Data_x, Table_2_Data_y,             linestyle='-',  color=color_2,    linewidth=LINE_WIDTH, label='2' )

        # plot analytical reference
        Sampling_x = np.linspace( 0.3*np.min(Table_1_Data_x), 3.0*np.max(Table_1_Data_x), 256 )

        ax.plot( Sampling_x,     Analytical_Ref(Sampling_x), linestyle='--', color=color_ref,  linewidth=LINE_WIDTH, label='Analytical' )

        # annotate the arrow and text
        #ax.annotate( '', xy=(annotated_arrow_1_x, annotated_arrow_1_y ), xytext=(annotated_text_1_x, annotated_text_1_y), va='bottom', ha='center', arrowprops=dict(arrowstyle='->', color=color_1, linewidth=LINE_WIDTH ) )
        #ax.text(annotated_text_1_x, annotated_text_1_y, annotated_text_1, va='center', ha='left', color=color_ref)

        #ax.annotate( '', xy=(annotated_arrow_2_x, annotated_arrow_2_y ), xytext=(annotated_text_2_x, annotated_text_2_y), va='bottom', ha='center', arrowprops=dict(arrowstyle='->', color=color_2, linewidth=LINE_WIDTH ) )
        #ax.text(annotated_text_2_x, annotated_text_2_y, annotated_text_2, va='center', ha='left', color=color_ref)

        # x,y scale
        ax.set_xscale( 'log' )
        ax.set_yscale( 'log' )

        # x,y labels
        if ( panel_idx//n_cols >= n_rows-1 ):
            ax.set_xlabel( axis_x_label )

        if ( panel_idx%n_cols == 0 ):
            ax.set_ylabel( axis_y_label )

        # x,y limit
        ax.set_xlim( axis_x_min, axis_x_max )
        ax.set_ylim( axis_y_min, axis_y_max )

        # legend
        if ( panel_idx == n_cols-1 ):
            ax.legend( loc='upper right' )

        # ticks and tick labels

        if ( panel_idx%n_cols == 0 ):
            ax.tick_params(axis='y', which='both', labelleft=True, labelright=False)
        else:
            ax.tick_params(axis='y', which='both', labelleft=False, labelright=False)

        if ( panel_idx//n_cols >= n_rows-1 ):
            ax.tick_params(axis='x', which='both', labelbottom=True, labeltop=False)
        else:
            ax.tick_params(axis='x', which='both', labelbottom=False, labeltop=False)

        minor_ticks = LogLocator( base=10.0, subs=np.arange(1.0, 10.0)*0.1, numticks=10 )
        ax.xaxis.get_ticklocs( minor=True )
        ax.yaxis.get_ticklocs( minor=True )
        ax.minorticks_on()
        ax.xaxis.set_minor_locator( minor_ticks )
        ax.xaxis.set_minor_formatter( NullFormatter() )
        ax.yaxis.set_minor_locator( minor_ticks )
        ax.yaxis.set_minor_formatter( NullFormatter() )
        ax.xaxis.set_ticks_position( 'both' )
        ax.yaxis.set_ticks_position( 'both' )
        ax.tick_params( which='both', direction='in' )

    # Output the figure to files
    plt.tight_layout(pad=0.2)
    fig.set_dpi(dpi)
    fig.set_size_inches(fig_size_x, fig_size_y)

    fig.savefig( 'fig_standard_paper_figure_1D.png', dpi=dpi )
    fig.savefig( 'fig_standard_paper_figure_1D.pdf', dpi=dpi )


if __name__ == '__main__':
    main()
