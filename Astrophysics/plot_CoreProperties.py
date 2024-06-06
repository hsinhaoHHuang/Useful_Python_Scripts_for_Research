import yt
import numpy as np
import matplotlib.pyplot as plt

ds0 = yt.load( '../Data_000000' )

# load the reference profiles
Step, Time, CoreDensity, CoreRadius, CoreMass, CoreRadiusCoreMass = np.loadtxt( '../Record__CoreProperties', skiprows=1, unpack=True )

# assign the units
Time                = ds0.arr( Time,               'code_time'             )
CoreDensity         = ds0.arr( CoreDensity,        'code_density'          )
CoreRadius          = ds0.arr( CoreRadius,         'code_length'           )
CoreMass            = ds0.arr( CoreMass,           'code_mass'             )
CoreRadiusCoreMass  = ds0.arr( CoreRadiusCoreMass, 'code_length*code_mass' )

# decide the units for plotting
#UNIT_L_PLOT = 'code_length'
#UNIT_D_PLOT = 'code_density'
UNIT_T_PLOT  = 'Myr'
UNIT_L_PLOT  = 'kpc'
UNIT_M_PLOT  = 'Msun'
UNIT_D_PLOT  = 'Msun/kpc**3'
UNIT_LM_PLOT = 'Msun*kpc'

# create the figure
fig = plt.figure()
ax1  = fig.add_subplot(221)
ax2  = fig.add_subplot(222)
ax3  = fig.add_subplot(223)
ax4  = fig.add_subplot(224)

# plot the profiles
ax1.plot( Time.in_units(UNIT_T_PLOT).d, CoreDensity.in_units(UNIT_D_PLOT).d,         color='C2', linewidth=0.3, label='Core Density' )
ax2.plot( Time.in_units(UNIT_T_PLOT).d, CoreRadius.in_units(UNIT_L_PLOT).d,          color='C3', linewidth=0.3, label='Core Radius'  )
ax3.plot( Time.in_units(UNIT_T_PLOT).d, CoreMass.in_units(UNIT_M_PLOT).d,            color='C0', linewidth=0.3, label='Core Mass'    )
ax4.plot( Time.in_units(UNIT_T_PLOT).d, CoreRadiusCoreMass.in_units(UNIT_LM_PLOT).d, color='C1', linewidth=0.3, label='CoreR*CoreM'  )

# set the limits and scales
# ax.set_xlim( xlim_min, xlim_max )
# ax.set_ylim( ylim_min, ylim_max )
# ax.set_xscale('log')
# ax.set_yscale('log')

# set the labels
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

ax1.set_ylabel( r'$\rho_c$'+r' ($M_{\odot}/{\rm kpc}^3$)' )
ax2.set_ylabel( r'$r_c$'+r' (${\rm kpc}$)' )
ax3.set_xlabel( r'$t$'+r' (${\rm Myr}$)'    )
ax3.set_ylabel( r'$M_c$'+r' ($M_{\odot}$)' )
ax4.set_xlabel( r'$t$'+r' (${\rm Myr}$)'    )
ax4.set_ylabel( r'$r_cM_c$'+r' ($M_{\odot}{\rm kpc}$)' )

# set the grid and ticks
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax1.tick_params( which='both',direction='in' )
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax2.tick_params( which='both',direction='in' )
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')
ax3.tick_params( which='both',direction='in' )
ax4.xaxis.set_ticks_position('both')
ax4.yaxis.set_ticks_position('both')
ax4.tick_params( which='both',direction='in' )

# save the figure
plt.tight_layout( pad=0.1, h_pad=0.1, w_pad=0.1 )
fig.savefig( 'fig_CoreProperties.png', dpi=150 )
fig.clear()
