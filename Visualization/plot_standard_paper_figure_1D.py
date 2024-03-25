#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np



####################################################
# Part 0. Setting
####################################################

####################################################
# Part 0.1 Genereal Setting
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
# Part 0.2 Setting for Plotting
####################################################
prof_color_FDM_16  = 'C3'       # color for FDM m_22=1.6 profile
prof_color_FDM_08  = 'C0'       # color for FDM m_22=1.6 profile
prof_color_CDM     = 'C2'       # color for FDM m_22=1.6 profile
prof_dens_min      = 3.9720e4   # upper limit for profile plot (in Msun/kpc**3)
prof_dens_max      = 2.3832e8   # lower limit for profile plot (in Msun/kpc**3)
prof_x_min         = 1.0e-1     # max r for profile plot (in kpc)
prof_x_max         = 3.0e1      # min r for profile plot (in kpc)

cm                 = 1/2.54     # centimeters in inches
fig_size_x         = 16*cm      # output figure size
fig_size_y         = 9*cm       # output figure size
fig                = plt.figure(1,(fig_size_x, fig_size_y), dpi=dpi)

###################################################################################################


####################################################
# Part 0.3 Read the data and parameters
####################################################
# path
Table_Path     = './'
Table_Filename = ''

# load the table
Table_Data_x, Table_Data_y = np.loadtxt( Table_Path+Table_Filename, unpack=True )

###################################################################################################


####################################################
# Part 0.4 Physical Constants and Variables
####################################################
##
Table_x_max = np.max(Table_Data_x)
Table_x_min = np.min(Table_Data_x)
Table_y_max = np.max(Table_Data_y)
Table_y_max = np.max(Table_Data_y)

###################################################################################################


####################################################
# Part 0.5 Define Functions
####################################################
def Analytical_Ref(x):
    y_analytical = np.sin(x)

    return y_analytical

###################################################################################################

####################################################
# Part 2. Profiles
####################################################

####################################################
# Part 2.1 Prepare data
####################################################
FDM_Units_L_plot = "kpc"
CDM_Units_L_plot = "kpccm"
####################################################################################################

Sampling_x = np.linspace(0.0, 2.0*np.pi, 100)

####################################################
# Part 2.2 Plot profiles
####################################################
ax_rect_left     = -0.01             # left boundary of the figure panels
ax_rect_bottom   =  0.02             # bottom boundary of the figure panels
ax_rect_width    =  0.97             # width of the figure panels
ax_rect_height   =  0.96             # height of the figure panels

ax = fig.add_axes( rect=[ax_rect_left, ax_rect_bottom, ax_rect_width, ax_rect_height] )

# plot data
ax.plot( Table_Data_x, Table_Data_y, '.-',     color=color_1, label='FDM ($m_{22}=1.6$)',  linewidth=LINE_WIDTH    )

# plot analytical reference
ax.plot( Sampling_x, Analytical_Ref(Sampling_x),  '--', color=color_ref,   label='Analytical Reference',  linewidth=LINE_WIDTH    )

# annotate the arrow and text
ax.annotate('', xy=(annotated_arrow_x, annotated_arrow_y ), xytext=(annotated_text_x, annotated_text_y), va='bottom', ha='center', arrowprops=dict(arrowstyle='->', color=prof_color_FDM_16,    linewidth=LINE_WIDTH) )
ax.text(annotated_text_x, annotated_text_y, annotated_text, va='center', ha='left', color=_color_ref)

# x,y scale
ax.set_xscale("log")
ax.set_yscale("log")

# x,y labels
ax.yaxis.set_label_position('left')
ax.set_xlabel(r'$\mathit{r}\ ({\rm kpc})}$')
ax.set_ylabel(r'$\mathit{\rho}\ ({\rm M}_{\odot}{\rm kpc}^{\rm -3})$')

# x,y limit
ax.set_xlim( 0.3*Table_x_min, 3.0*Table_x_max )
ax.set_ylim( 0.3*Table_y_min, 3.0*Table_y_max )

# legend
ax.legend( loc='upper right' )

# ticks and tick labels
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(which='both',direction='in')
ax.tick_params(axis='y', which='both', labelleft=False, labelright=True)

####################################################################################################
####################################################################################################


####################################################
# Part 3. Output files
####################################################
fig.set_dpi(dpi)
fig.set_size_inches(fig_size_x, fig_size_y)

fig.savefig("fig.png", dpi=dpi)
fig.savefig("fig.pdf", dpi=dpi)
