import numpy as np
import matplotlib.pyplot as plt


###############################################################################
# General setting
dpi        = 150
FONT_SIZE  = 14.0
LINE_WIDTH = 1.0
MARK_SIZE  = 4.0

plt.rcParams['font.size']         = FONT_SIZE
plt.rcParams['figure.titlesize']  = FONT_SIZE
plt.rcParams['figure.dpi']        = dpi

plt.rcParams['axes.titlesize']    = FONT_SIZE
plt.rcParams['axes.labelsize']    = FONT_SIZE
plt.rcParams['axes.labelpad']     = FONT_SIZE/8.0*0.05
plt.rcParams['axes.linewidth']    = LINE_WIDTH

plt.rcParams['legend.fontsize']   = FONT_SIZE
plt.rcParams['lines.linewidth']   = LINE_WIDTH
plt.rcParams["legend.labelspacing"] = 0.02
plt.rcParams["legend.borderpad"]  = 0.1

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

###############################################################################
# Load the table
Table_filename = "Data_000000_halo#1_condensation_profile_au"

with open(Table_filename, 'r') as f:
    for line in f:
        if line.startswith('#'):
            header = line
        else:
            break #stop when there are no more #

Table_header = header[1:].strip().split()
Table_data   = np.genfromtxt( Table_filename, delimiter=None, comments='#', names=Table_header, dtype=None, encoding=None )

###############################################################################

VelocityDispersion = np.fromfile("velocity_dispersion.bin")
VelocityDispersion_r     = VelocityDispersion[:VelocityDispersion.shape[0]//2]
VelocityDispersion_sigma = VelocityDispersion[VelocityDispersion.shape[0]//2:]

###############################################################################
Mpc_h_to_kpc    = 1.48541685e+03
m_to_km         = 1e-3
km_to_m         = 1e+3
codemass_to_kg  = 2.57725326e+41

Newton_G        = 4.30078846e-06
ELBDM_Eta       = 7.74803816e+03/100.0/Mpc_h_to_kpc # in (s/km)(1/kpc)
ELBDM_Mass      = 1.0e-22                           # in eV/c^2

###############################################################################
def Diameter_granule( sigma ): # in kpc (sigma in km/s)
   return ( (36.0*np.pi)**(1.0/6.0) )/(ELBDM_Eta*sigma)
   #return 0.35*(2.0*np.pi)/(ELBDM_Eta*sigma)

def kbTemperature_granule( sigma ): # in eV (sigma in km/s)
   return ELBDM_Mass*( (sigma/3.0e5)**2 )

###############################################################################
# Plot the data
cm          = 1/2.54     # centimeters in inches
fig_size_x  = 16*cm      # output figure size
fig_size_y  = 18*cm      # output figure size
fig, ax     = plt.subplots(3, 1, figsize=(fig_size_x, fig_size_y), dpi=dpi)

Orbital_v = np.sqrt(Newton_G*Table_data["Condensation_Env_M"]/Table_data["Condensation_Env_R"])

ax[0].plot( Table_data["Condensation_Env_R"],  Table_data["Condensation_v"]/np.sqrt(3),        '-o',  color='C0',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{\frac{1}{3}\langle v^2\rangle_{\rm enc}}$"    )
ax[0].plot( Table_data["Condensation_Env_R"],  Table_data["Condensation_Shell_v"]/np.sqrt(3),  '-o',  color='C1',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{\frac{1}{3}\langle v^2\rangle_{\rm shell}}$"  )
#ax[0].plot( Table_data["Condensation_Env_R"],  Orbital_v/np.sqrt(2),                           '-o',  color='C2',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{GM_{\rm enc}/r}/\sqrt{2}$"                    )
ax[0].plot( VelocityDispersion_r*Mpc_h_to_kpc, VelocityDispersion_sigma*m_to_km*np.sqrt(2),    '--',  color='C3',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r'$\sigma_{\rm Jeans}$'                                )

# Set the limit
ax[0].set_xlim( Table_data["Condensation_Env_R"][0]-1 , Table_data["Condensation_Env_R"][-1]+1 )
ax[0].set_ylim( 10.0, 34.0 )

# Set the label
ax[0].set_xlabel( r"$r\ ({\rm kpc})$" )
ax[0].set_ylabel( r"$\sigma\ ({\rm km/s})$" )

# Set the legend
ax[0].legend(loc = 'upper right')

# Set grids and ticks
ax[0].grid()
ax[0].xaxis.set_ticks_position('both')
ax[0].yaxis.set_ticks_position('both')
ax[0].minorticks_on()
ax[0].tick_params(which='both',direction='in')

###########################################################################

ax[1].plot( Table_data["Condensation_Env_R"],  Diameter_granule( Table_data["Condensation_v"]/np.sqrt(3) ),        '-o',  color='C0',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{\frac{1}{3}\langle v^2\rangle_{\rm enc}}$"   )
ax[1].plot( Table_data["Condensation_Env_R"],  Diameter_granule( Table_data["Condensation_Shell_v"]/np.sqrt(3) ),  '-o',  color='C1',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{\frac{1}{3}\langle v^2\rangle_{\rm shell}}$" )
#ax[1].plot( Table_data["Condensation_Env_R"],  Diameter_granule( Orbital_v/np.sqrt(2) ),                           '-o',  color='C2',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{GM_{\rm enc}/r}/\sqrt{2}$"                   )
ax[1].plot( VelocityDispersion_r*Mpc_h_to_kpc, Diameter_granule( VelocityDispersion_sigma*m_to_km*np.sqrt(2) ),    '--',  color='C3',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r'$\sigma_{\rm Jeans}$'                               )

# Set the limit
ax[1].set_xlim( Table_data["Condensation_Env_R"][0]-1 , Table_data["Condensation_Env_R"][-1]+1 )
ax[1].set_ylim( 1.0, 4.2 )
#ax[1].set_yscale("log")

# Set the label
ax[1].set_xlabel( r"$r\ ({\rm kpc})$" )
ax[1].set_ylabel( r"$d_{\rm granule}\ ({\rm kpc})$" )

# Set grids and ticks
ax[1].grid()
ax[1].xaxis.set_ticks_position('both')
ax[1].yaxis.set_ticks_position('both')
ax[1].minorticks_on()
ax[1].tick_params(which='both',direction='in')

###########################################################################

ax[2].plot( Table_data["Condensation_Env_R"],  kbTemperature_granule( Table_data["Condensation_v"]/np.sqrt(3) ),        '-o',  color='C0',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{\frac{1}{3}\langle v^2\rangle_{\rm enc}}$"   )
ax[2].plot( Table_data["Condensation_Env_R"],  kbTemperature_granule( Table_data["Condensation_Shell_v"]/np.sqrt(3) ),  '-o',  color='C1',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{\frac{1}{3}\langle v^2\rangle_{\rm shell}}$" )
#ax[2].plot( Table_data["Condensation_Env_R"],  kbTemperature_granule( Orbital_v/np.sqrt(2) ),                           '-o',  color='C2',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r"$\sqrt{GM_{\rm enc}/r}/\sqrt{2}$"                   )
ax[2].plot( VelocityDispersion_r*Mpc_h_to_kpc, kbTemperature_granule( VelocityDispersion_sigma*m_to_km*np.sqrt(2) ),    '--',  color='C3',  linewidth=LINE_WIDTH,  markersize=MARK_SIZE,  label=r'$\sigma_{\rm Jeans}$'                               )

# Set the limit
ax[2].set_xlim( Table_data["Condensation_Env_R"][0]-1 , Table_data["Condensation_Env_R"][-1]+1 )
ax[2].set_ylim( 1.0e-31, 1.2e-30 )

# Set the label
ax[2].set_xlabel( r"$r\ ({\rm kpc})$" )
ax[2].set_ylabel( r"$k_bT_{\rm granule}\ ({\rm eV})$" )

# Set grids and ticks
ax[2].grid()
ax[2].xaxis.set_ticks_position('both')
ax[2].yaxis.set_ticks_position('both')
ax[2].minorticks_on()
ax[2].tick_params(which='both',direction='in')

###########################################################################

# Save the figure
fig.set_dpi(dpi)
fig.set_size_inches( fig_size_x, fig_size_y )

plt.tight_layout( pad=0.1, h_pad=0.1, w_pad=0.1 )
fig.savefig("fig_GranuleSize.png", dpi=dpi)
fig.clear()
