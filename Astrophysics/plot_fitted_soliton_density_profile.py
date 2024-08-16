import argparse
import sys
import yt
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Profiles' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start    = args.idx_start
idx_end      = args.idx_end
didx         = args.didx
prefix       = args.prefix

r_sphere     = (100.0, 'kpc')
dpi          = 150
nbin         = 32

Ref_DensProf_filename = '../HaloDensityProfile'

yt.enable_parallelism()
ts = yt.DatasetSeries([ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ])

###################################################################################################################
# code units (in cgs)
UNIT_L     = 4.4366320345075490e+24
UNIT_D     = 2.5758579724476994e-30

# length units
Const_cm   = 1.0;
Const_pc   = 3.08567758149e18;
Const_kpc  = 1.0e3*Const_pc;

# mass units
Const_g    = 1.0;
Const_kg   = 1.0e3*Const_g;
Const_Msun = 1.9885e33;

def Soliton_rho_c_from_r_c(m22, r_c):
    return 1.945e7/( m22**2 * (r_c*UNIT_L/Const_kpc)**4 )*(Const_Msun/Const_kpc**3)/UNIT_D

def Soliton_fitting_analytical_dens(r, m22, r_c):
    rho_c = Soliton_rho_c_from_r_c(m22, r_c)
    x = r / r_c
    return rho_c*( 1.0+9.06e-2*(x**2) )**(-8)
###################################################################################################################

for ds in ts.piter():

   fields_list = ['density']

   if ds.parameters["Particle"] == 1:
      if ds.parameters["Par_NPar"] > 0:
         fields_list.append('particle_density_on_grid')

   for field in fields_list:

      # set the center as the position of the gravitational potential minimum
      center_pos  = ds.all_data().quantities.min_location(('gamer','Pote'))[1:]
      sp          = ds.sphere( center_pos, r_sphere )

      # ProfilePlot
      prof = yt.ProfilePlot( sp, 'radius', field, weight_field='cell_volume', n_bins=nbin )
      prof.set_unit( 'radius', 'kpc' )
      prof.set_xlim(        1.0e-1, 1.0e+2 )
      prof.set_ylim( field, 1.0e+2, 2.0e+7 )
      prof.set_unit( field, 'Msun/kpc**3'  )
      prof.annotate_title( 't = %13.7e Gyr'%(ds.current_time.in_units('Gyr')) )

      prof.save( 'fig_%s.png'%(ds), mpl_kwargs={'dpi':dpi} )

      # create the arrays of profile
      prof_dens = yt.create_profile( sp, 'radius', fields=field,
                                     weight_field='cell_volume', n_bins=nbin,
                                     units={'radius':'code_length', field:'code_density'} )

      # save the profile to text file
      np.savetxt( '%s_%s_profile'%(ds,field),
                  np.column_stack( (prof_dens.x.in_units('code_length').d, prof_dens[field].in_units('code_density').d) ),
                  fmt='          %9.8e',
                  header='       r (code_length)   density (code_density)' )

      # Find the core radius by curve fit

      m22             = 1.0
      peak_dens       = prof_dens[field].in_units('code_density').d[0]
      central_region  = (prof_dens[field].in_units('code_density').d > 0.1*peak_dens)
      popt, pcov      = scipy.optimize.curve_fit( lambda x, r_c: Soliton_fitting_analytical_dens( x, m22, r_c ),
                                                  prof_dens.x.in_units('code_length').d[central_region],
                                                  prof_dens[field].in_units('code_density').d[central_region],
                                                  p0=prof_dens.x.in_units('code_length').d[central_region][-1],
                                                  bounds=(0, np.inf) )
      core_radius           = ds.arr( popt[0], 'code_length' )
      core_density          = ds.arr( Soliton_rho_c_from_r_c(m22, core_radius.in_units('code_length').d ), 'code_density' )
      Soliton_DensProf_r    = ds.arr(  prof_dens.x.in_units('code_length').d, 'code_length' )
      Soliton_DensProf_dens = ds.arr( Soliton_fitting_analytical_dens( Soliton_DensProf_r.in_units('code_length').d, m22, core_radius.in_units('code_length').d ), 'code_density' )

      # load the reference profiles
      Ref_DensProf_r, Ref_DensProf_dens = np.loadtxt( Ref_DensProf_filename,                 skiprows=1, unpack=True )
      Ini_DensProf_r, Ini_DensProf_dens = np.loadtxt( '%s_%s_profile'%('Data_000000',field), skiprows=1, unpack=True )

      # assign the units
      Ref_DensProf_r    = ds.arr( Ref_DensProf_r,    'code_length'  )
      Ref_DensProf_dens = ds.arr( Ref_DensProf_dens, 'code_density' )
      Ini_DensProf_r    = ds.arr( Ini_DensProf_r,    'code_length'  )
      Ini_DensProf_dens = ds.arr( Ini_DensProf_dens, 'code_density' )

      # decide the units for plotting
      #UNIT_L_PLOT = 'code_length'
      #UNIT_D_PLOT = 'code_density'
      UNIT_L_PLOT = 'kpc'
      UNIT_D_PLOT = 'Msun/kpc**3'

      # create the figure
      fig = plt.figure()
      ax  = fig.add_subplot(111)

      # plot the profiles
      ax.plot( Ref_DensProf_r.in_units(UNIT_L_PLOT).d, Ref_DensProf_dens.in_units(UNIT_D_PLOT).d, label='Reference'   )
      ax.plot( Ini_DensProf_r.in_units(UNIT_L_PLOT).d, Ini_DensProf_dens.in_units(UNIT_D_PLOT).d, label='Data_000000' )
      ax.plot( prof_dens.x.in_units(UNIT_L_PLOT).d,    prof_dens[field].in_units(UNIT_D_PLOT).d,  label='%s'%ds       )
      ax.plot(    Soliton_DensProf_r.in_units(UNIT_L_PLOT).d, Soliton_DensProf_dens.in_units(UNIT_D_PLOT).d, linestyle='--', c='red',  label='Soliton'                                                                           )
      ax.axvline( core_radius.in_units(UNIT_L_PLOT).d,                                                       linestyle='-.', c='pink', label=r'$r_c$=%4.2e %s'%(core_radius.in_units(UNIT_L_PLOT).d, UNIT_L_PLOT ),     zorder=1 )
      ax.axhline( core_density.in_units(UNIT_D_PLOT).d,                                                      linestyle=':',  c='pink', label=r'$\rho_c$=%4.2e %s'%(core_density.in_units(UNIT_D_PLOT).d, UNIT_D_PLOT ), zorder=1 )

      # set the limits and scales
      ax.set_xlim( 0.5*np.min(Ref_DensProf_r.in_units(UNIT_L_PLOT).d),     2.0*np.max(Ref_DensProf_r.in_units(UNIT_L_PLOT).d)    )
      ax.set_ylim( 0.5*np.min(Ref_DensProf_dens.in_units(UNIT_D_PLOT).d), 20.0*np.max(Ref_DensProf_dens.in_units(UNIT_D_PLOT).d) )
      ax.set_xscale('log')
      ax.set_yscale('log')

      # set the labels
      ax.legend()
      ax.set_xlabel( r'$r$'+' (%s)'%UNIT_L_PLOT    )
      ax.set_ylabel( r'$\rho$'+' (%s)'%UNIT_D_PLOT )
      ax.set_title( '$t$ = %7.6e Gyr'%ds.current_time.in_units('Gyr') )

      # set the grid and ticks
      ax.grid()
      ax.xaxis.set_ticks_position('both')
      ax.yaxis.set_ticks_position('both')
      ax.tick_params( which='both',direction='in' )

      # save the figure
      plt.tight_layout( pad=0.1, h_pad=0.1, w_pad=0.1 )
      fig.savefig( 'fig_%s_Profile_%s.png'%(ds, field), dpi=dpi )
      fig.clear()
