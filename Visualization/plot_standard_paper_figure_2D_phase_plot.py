#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import yt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.ticker import LogLocator, NullFormatter


###################################################################################################
# Genereal Setting
###################################################################################################
dpi        = 600
FONT_SIZE  = 8.0
LINE_WIDTH = 1.0

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


###################################################################################################
# Setting for Plotting
###################################################################################################
n_rows             = 2                     # number of rows in the figure panels
n_cols             = 3                     # number of columns in the figure panels

# fields to plot the phase plots
field_z            = ('gas','cell_mass')   # target field
field_x            = ('gas','density')     # field at the x-axis
field_y            = ('gas','temperature') # field at the y-axis

Plotting_UNIT_L    = 'kpc'                 # length unit in the plot
Plotting_UNIT_D    = 'Msun/kpc**3'         # density unit in the plot
Plotting_UNIT_M    = 'Msun'                # density unit in the plot
Plotting_UNIT_T    = 'Myr'                 # time unit in the plot
Plotting_UNIT_Temp = 'K'                   # temperature unit in the plot

UNIT_D_todisplay   = r'$\mathit{\rho}\ ({\rm M}_{\odot}{\rm kpc}^{\rm -3})$'
UNIT_M_todisplay   = r'$\mathit{M}\ ({\rm M}_{\odot})$'
UNIT_Temp_todisplay= r'$\mathit{T}\ ({\rm K})$'

sphere_radius      = -1.0                  # radius of sphere for the phase plots (in plotting unit) (<=0.0=one-half box width)

colormap           = 'viridis'             # ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
color_lim_min      = 1.e0                  # color bar lower limit for z-field (in plotting unit)
color_lim_max      = 1.e7                  # color bar upper limit for z-field (in plotting unit)

x_lim_min          = 9.e1                  # lower limit for the x-field
x_lim_max          = 2.e5                  # upper limit for the x-field
y_lim_min          = 9.e2                  # lower limit for the y-field
y_lim_max          = 2.e6                  # upper limit for the y-field

grid_rect_left     = 0.05                  # left boundary of the figure panels
grid_rect_bottom   = 0.08                  # bottom boundary of the figure panels
grid_rect_width    = 0.90                  # width of the figure panels
grid_rect_height   = 0.90                  # height of the figure panels

cm                 = 1/2.54                # centimeters in inches
fig_size_x         = 16*cm                 # output figure size in the horizontal direction
fig_size_y         = (fig_size_x-1*cm)/n_cols*n_rows
fig                = plt.figure( 1,(fig_size_x, fig_size_y), dpi=dpi )
###################################################################################################


###################################################################################################
# Read the data and parameters
###################################################################################################
# path
Data_path    = np.empty( n_rows*n_cols, dtype=object )
Data_path[0] = './Data_000001'
Data_path[1] = './Data_000002'
Data_path[2] = './Data_000003'
Data_path[3] = './Data_000004'
Data_path[4] = './Data_000005'
Data_path[5] = './Data_000006'
###################################################################################################



def main() -> None:

    # Create the AxesGrid
    grid = AxesGrid( fig, rect=[grid_rect_left, grid_rect_bottom, grid_rect_width, grid_rect_height],
                     nrows_ncols=(n_rows, n_cols),
                     axes_pad=(0.05,0.05), label_mode="L", share_all=True,
                     cbar_location="right", cbar_mode="single", cbar_size="2%", cbar_pad="1%", aspect=False )

    # Loop for each panel
    for panel_idx in range(0, n_rows*n_cols, 1):

        ###################################################################################################
        # Load the dataset
        ds = yt.load( Data_path[panel_idx] )

        ###################################################################################################
        # Region to plot
        Data_Center_x         = ds.all_data().quantities.max_location('density')[1].in_units('code_length').d if sphere_radius > 0.0 else ds.domain_center[0].in_units('code_length').d
        Data_Center_y         = ds.all_data().quantities.max_location('density')[2].in_units('code_length').d if sphere_radius > 0.0 else ds.domain_center[1].in_units('code_length').d
        Data_Center_z         = ds.all_data().quantities.max_location('density')[3].in_units('code_length').d if sphere_radius > 0.0 else ds.domain_center[2].in_units('code_length').d
        Data_Center           = ds.arr( [Data_Center_x, Data_Center_y, Data_Center_z], 'code_length' )
        Data_sphere_radius    = ds.quan( sphere_radius, Plotting_UNIT_L ) if sphere_radius > 0.0 else 0.5*ds.domain_width[0]
        Data_region           = ds.sphere( Data_Center, Data_sphere_radius )

        ###################################################################################################
        # Phase plot
        phas = yt.PhasePlot( Data_region, field_x, field_y, field_z,
                                          weight_field=None, x_bins=300, y_bins=300 )

        # Set the font
        phas.set_font( { "size":FONT_SIZE, "math_fontfamily":'custom' } )

        # Set the field unit
        phas.set_unit( field_x, Plotting_UNIT_D    )
        phas.set_unit( field_y, Plotting_UNIT_Temp )
        phas.set_unit( field_z, Plotting_UNIT_M    )

        # Set the colorbar limits
        phas.set_xlim(              x_lim_min,     x_lim_max )
        phas.set_ylim(              y_lim_min,     y_lim_max )
        phas.set_zlim( field_z, color_lim_min, color_lim_max )

        # Set the colormap
        phas.set_cmap( field_z, colormap )

        if ( panel_idx%n_cols == 0  and  panel_idx//n_cols >= n_rows-1 ):
            # Set the x label
            phas.set_xlabel( UNIT_D_todisplay )

            # Set the y label
            phas.set_ylabel( UNIT_Temp_todisplay )

        else:
            phas.set_xlabel( ' ' )
            phas.set_ylabel( ' ' )

        # Set the colorbar label
        phas.set_colorbar_label( field_z, UNIT_M_todisplay )

        # Annotate the time
        phas.annotate_text( xpos=x_lim_min*10**(0.60*(np.log10(x_lim_max)-np.log10(x_lim_min))),
                            ypos=y_lim_min*10**(0.90*(np.log10(y_lim_max)-np.log10(y_lim_min))),
                            text='$t$ = {:.1f} {:s}'.format( ds.current_time.in_units(Plotting_UNIT_T).d, Plotting_UNIT_T ),
                            size=FONT_SIZE, color="black" )


        ###################################################################################################
        # Set the plots to the grid
        plot        = phas.plots[field_z]
        plot.figure = fig
        plot.axes   = grid[panel_idx].axes
        plot.cax    = grid.cbar_axes[0]

        # actually render the phase plot
        phas.render()

        # Set the ticks
        minor_ticks = LogLocator( base=10.0, subs=np.arange(1.0, 10.0)*0.1, numticks=10 )
        grid[panel_idx].axes.xaxis.get_ticklocs( minor=True )
        grid[panel_idx].axes.yaxis.get_ticklocs( minor=True )
        grid[panel_idx].axes.minorticks_on()
        grid[panel_idx].axes.xaxis.set_minor_locator( minor_ticks )
        grid[panel_idx].axes.xaxis.set_minor_formatter( NullFormatter() )
        grid[panel_idx].axes.yaxis.set_minor_locator( minor_ticks )
        grid[panel_idx].axes.yaxis.set_minor_formatter( NullFormatter() )
        grid[panel_idx].axes.tick_params( which='both', left=True, right=True, bottom=True, top=True )
        grid[panel_idx].tick_params( which='both', direction='in' )

        # only the L-shape
        if ( panel_idx%n_cols == 0 ):
            grid[panel_idx].axes.tick_params( labelleft=True, labelright=False )
        else:
            grid[panel_idx].axes.tick_params( labelleft=False, labelright=False )

        if ( panel_idx//n_cols >= n_rows-1 ):
            grid[panel_idx].axes.tick_params( labelbottom=True, labeltop=False )
        else:
            grid[panel_idx].axes.tick_params( labelbottom=False, labeltop=False )

        # Annotate the text for title
        grid[panel_idx].axes.text( 0.05, 0.90, 'Data %d'%(panel_idx), transform=grid[panel_idx].axes.transAxes,
                                   size=1.25*FONT_SIZE, color="black", weight="normal" )

        ###################################################################################################

    # Output the figure to files
    fig.set_dpi( dpi )
    fig.set_size_inches( fig_size_x, fig_size_y )

    fig.savefig( "fig_standard_paper_figure_2D_phase_plot.png", dpi=dpi )
    fig.savefig( "fig_standard_paper_figure_2D_phase_plot.pdf", dpi=dpi )


if __name__ == '__main__':
    main()
