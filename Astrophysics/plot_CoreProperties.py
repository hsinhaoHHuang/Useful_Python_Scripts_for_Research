import yt
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

ds0 = yt.load( '../Data_000000' )

def X_IntersectionInLogSpace( X_1, X_2, F_1, F_2, G_1, G_2  ):
    # Linear interpolation in log-space
    x_1 = np.log10( X_1 )
    x_2 = np.log10( X_2 )
    f_1 = np.log10( F_1 )
    f_2 = np.log10( F_2 )
    g_1 = np.log10( G_1 )
    g_2 = np.log10( G_2 )

    P_X = 10**( (x_1*f_2 - x_2*f_1 - x_1*g_2 + x_2*g_1)/(g_1 - g_2 - f_1+ f_2) )
    return P_X

# load the reference profiles
Step, Time, CoreDensity, CoreRadius, CoreMass, CoreRadiusCoreMass = np.loadtxt( '../Record__CoreProperties', skiprows=1, unpack=True )

# decide the units for plotting
#UNIT_L_PLOT = 'code_length'
#UNIT_D_PLOT = 'code_density'
UNIT_T_PLOT  = 'Myr'
UNIT_L_PLOT  = 'kpc'
UNIT_M_PLOT  = 'Msun'
UNIT_D_PLOT  = 'Msun/kpc**3'
UNIT_LM_PLOT = 'Msun*kpc'

# assign the units
Time                = ds0.arr( Time,               'code_time'             ).in_units(UNIT_T_PLOT).d
CoreDensity         = ds0.arr( CoreDensity,        'code_density'          ).in_units(UNIT_D_PLOT).d
CoreRadius          = ds0.arr( CoreRadius,         'code_length'           ).in_units(UNIT_L_PLOT).d
CoreMass            = ds0.arr( CoreMass,           'code_mass'             ).in_units(UNIT_M_PLOT).d
CoreRadiusCoreMass  = ds0.arr( CoreRadiusCoreMass, 'code_length*code_mass' ).in_units(UNIT_LM_PLOT).d

# number of points for interpolation
N_interp                       = 100*(Time.size - 1)                   # 100 times the original data
dt_interp                      = (Time[-1]-Time[0])/(N_interp)
Time_interp                    = np.arange( Time[0], Time[-1], dt_interp )

# interpolation to make sure it is evenly spaced
CoreDensity_interp             = np.interp( Time_interp, Time, CoreDensity        )
CoreRadius_interp              = np.interp( Time_interp, Time, CoreRadius         )
CoreMass_interp                = np.interp( Time_interp, Time, CoreMass           )
CoreRadiusCoreMass_interp      = np.interp( Time_interp, Time, CoreRadiusCoreMass )

# Gaussian sigma for smoothing
Soliton_tau_00                 = 39.9*(np.min(CoreDensity)/1.0e9)**(-0.5) # Soliton wavefunction oscillation period
Gaussian_sigma_to_tau00_ratio  = 8
Gaussian_sigma                 = Gaussian_sigma_to_tau00_ratio*Soliton_tau_00
Gaussian_truncate              = 5.0           # in sigma

# Gaussian smoothing
CoreDensity_smooth_filt        = gaussian_filter1d( CoreDensity_interp,        sigma=Gaussian_sigma/dt_interp, mode='nearest', truncate=Gaussian_truncate )
CoreRadius_smooth_filt         = gaussian_filter1d( CoreRadius_interp,         sigma=Gaussian_sigma/dt_interp, mode='nearest', truncate=Gaussian_truncate )
CoreMass_smooth_filt           = gaussian_filter1d( CoreMass_interp,           sigma=Gaussian_sigma/dt_interp, mode='nearest', truncate=Gaussian_truncate )
CoreRadiusCoreMass_smooth_filt = gaussian_filter1d( CoreRadiusCoreMass_interp, sigma=Gaussian_sigma/dt_interp, mode='nearest', truncate=Gaussian_truncate )

# define the threshold value for soliton formation
CoreDensity_threshold_ratio        = 2.0
CoreDensity_threshold_value        = CoreDensity_threshold_ratio*CoreDensity[0]
CoreDensity_threshold_array        = np.full( np.shape(Time_interp), CoreDensity_threshold_value )

CoreRadiusCoreMass_threshold_ratio = 0.9
CoreRadiusCoreMass_threshold_value = CoreRadiusCoreMass_threshold_ratio*CoreRadiusCoreMass[0]
CoreRadiusCoreMass_threshold_array = np.full( np.shape(Time_interp), CoreRadiusCoreMass_threshold_value )

# find the intersection of smooth curve and the threshold value
# peak density intersection
CoreDensity_intersection_idx          = np.argwhere( np.diff(np.sign(CoreDensity_threshold_array - CoreDensity_smooth_filt)) ).flatten()[0]
try:
   CoreDensity_intersection_Time      = X_IntersectionInLogSpace( Time_interp[CoreDensity_intersection_idx],
                                                                  Time_interp[CoreDensity_intersection_idx+1],
                                                                  CoreDensity_threshold_array[CoreDensity_intersection_idx],
                                                                  CoreDensity_threshold_array[CoreDensity_intersection_idx+1],
                                                                  CoreDensity_smooth_filt[CoreDensity_intersection_idx],
                                                                  CoreDensity_smooth_filt[CoreDensity_intersection_idx+1] )
except Exception as e:
   print(e)
   CoreDensity_intersection_Time      = np.nan
   pass

# rcmc intersection
CoreRadiusCoreMass_intersection_idx         = np.argwhere( np.diff(np.sign(CoreRadiusCoreMass_threshold_array - CoreRadiusCoreMass_smooth_filt)) ).flatten()[0]
try:
   CoreRadiusCoreMass_intersection_Time     = X_IntersectionInLogSpace( Time_interp[CoreRadiusCoreMass_intersection_idx],
                                                                        Time_interp[CoreRadiusCoreMass_intersection_idx+1],
                                                                        CoreRadiusCoreMass_threshold_array[CoreRadiusCoreMass_intersection_idx],
                                                                        CoreRadiusCoreMass_threshold_array[CoreRadiusCoreMass_intersection_idx+1],
                                                                        CoreRadiusCoreMass_smooth_filt[CoreRadiusCoreMass_intersection_idx],
                                                                        CoreRadiusCoreMass_smooth_filt[CoreRadiusCoreMass_intersection_idx+1] )
except Exception as e:
   print(e)
   CoreRadiusCoreMass_intersection_Time     = np.nan
   pass

###################################################################################################

# create the figure
fig = plt.figure()
ax1  = fig.add_subplot(221)
ax2  = fig.add_subplot(222)
ax3  = fig.add_subplot(223)
ax4  = fig.add_subplot(224)

# plot the profiles
ax1.plot( Time, CoreDensity,         color='C2', linewidth=0.3, label='Core Density' )
ax2.plot( Time, CoreRadius,          color='C3', linewidth=0.3, label='Core Radius'  )
ax3.plot( Time, CoreMass,            color='C0', linewidth=0.3, label='Core Mass'    )
ax4.plot( Time, CoreRadiusCoreMass,  color='C1', linewidth=0.3, label='CoreR*CoreM'  )

ax1.plot( Time_interp, CoreDensity_smooth_filt,     '+-',  color='C3',   markersize=4,  label="smooth, nearest, $\sigma$=%.2e Myr"%Gaussian_sigma  )
ax1.plot( Time_interp, CoreDensity_threshold_array, '--',  color='grey', markersize=4,  label="threshold value (%.2f)"%CoreDensity_threshold_ratio )

if not np.isnan(CoreDensity_intersection_Time):
    ax1.plot( [CoreDensity_intersection_Time, CoreDensity_intersection_Time], [np.min(CoreDensity), np.max(CoreDensity)], '-.',  color='r', markersize=4,  label="_Time, cross the threshold" )
    ax1.annotate( '%7.6e Myr'%CoreDensity_intersection_Time, xy=(CoreDensity_intersection_Time, CoreDensity_threshold_value ),
                   xytext=(1.5*CoreDensity_intersection_Time, CoreDensity[0] ),
                   va='bottom', ha='center', arrowprops=dict(arrowstyle='->', color='k', linewidth=1.5) )

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
