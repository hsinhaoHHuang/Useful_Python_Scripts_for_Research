import yt
import numpy as np


def Add_profile_normalized_field( ds, field_name, center, r_max ):
    sp             = ds.sphere( center, r_max )
    profile        = yt.create_profile( sp, 'radius', fields=field_name, weight_field='cell_volume', n_bins=32 )
    profile_radius = profile.x
    profile_value  = profile[field_name]

    def _profile_normalized_field( field, data ):
        r = ( (data[('gamer','x')]-center[0])**2 + (data[('gamer','y')]-center[1])**2 + (data[('gamer','z')]-center[2])**2 )**0.5
        return data[field_name]/( np.interp( r/profile_radius.uq, profile_radius/profile_radius.uq, profile_value/profile_value.uq )*profile_value.uq )

    ds.add_field(       ('gamer', 'profile_normalized_'+field_name[1] ),
                  function      = _profile_normalized_field,
                  display_name  = 'Profile Normalized '+' '.join(field_name[1].split('_')),
                  units         = 'dimensionless',
                  sampling_type = 'cell' )


def Add_reciprocal_field( ds, field_name ):

    def _reciprocal_field( field, data ):
        return data[field_name]**(-1)

    ds.add_field(       ('gamer', 'reciprocal_'+field_name[1] ),
                  function      = _reciprocal_field,
                  display_name  = 'Reciprocal '+' '.join(field_name[1].split('_')),
                  sampling_type = 'cell' )
