#!/usr/bin/env python3

import yt
import matplotlib.pyplot as plt

Filename_1 = ''
Filename_2 = ''

def max_density_position( dataset ):

    center_pos  = dataset.all_data().quantities.max_location('density')[1:]

    return center_pos

def half_box_size_x( dataset ):

    half_Lx = 0.5*dataset.domain_width[0]

    return half_Lx

def density_profile( sphere, number_of_bins: int ):

    dens_prof = yt.create_profile( sphere, 'radius', fields='density', weight_field='cell_volume', n_bins=number_of_bins )

    return dens_prof

def arrays_of_density_profile( density_profile, unit_radius: str, unit_density: str ) -> tuple[float, float]:

    array_radius   = density_profile.x.in_units(unit_radius).d
    array_density  = density_profile['density'].in_units(unit_density).d

    return array_radius, array_density

def plot_two_density_profiles_to_figure( radius_1: float, density_1: float, radius_2: float, density_2:float ) -> None:

    # create the figure
    fig = plt.figure()
    ax  = fig.add_subplot(111)

    # plot the profiles
    ax.plot( radius_1, density_1, label='Density profile 1' )
    ax.plot( radius_2, density_2, label='Density profile 2' )

    # set scales
    ax.set_xscale('log')
    ax.set_yscale('log')

    # set lables
    ax.legend()
    ax.set_xlabel( r'$r\ ({\rm kpc})$'              )
    ax.set_ylabel( r'$\rho\ (M_{\odot}/{\rm kpc}^3)$' )

    # save the figure
    plt.tight_layout()
    fig.savefig( 'fig_TwoDensityProfiles.png' )

def main() -> None:
    # load the datasets
    ds_1 = yt.load( Filename_1 )
    ds_2 = yt.load( Filename_2 )

    # create the spheres, using peak density as center and half box size as radius
    sphere_1 = ds_1.sphere( max_density_position(ds_1), half_box_size_x(ds_1) )
    sphere_2 = ds_2.sphere( max_density_position(ds_2), half_box_size_x(ds_2) )

    # create the density profiles, using 32 bins, with 'kpc' as length unit, and 'Msun/kpc**3' as density unit
    dens_prof_radius_1, dens_prof_density_1 = arrays_of_density_profile( density_profile( sphere_1, 32 ), 'kpc', 'Msun/kpc**3' )
    dens_prof_radius_2, dens_prof_density_2 = arrays_of_density_profile( density_profile( sphere_2, 32 ), 'kpc', 'Msun/kpc**3' )

    # plot the density profiles
    plot_two_density_profiles_to_figure( dens_prof_radius_1, dens_prof_density_1, dens_prof_radius_2, dens_prof_density_2 )

if __name__ == '__main__':
    main()
