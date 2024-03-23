#!/usr/bin/env python3

import numpy as np

# Set the print formatter for float arrays
float_formatter = '{: 6.1f}'.format
np.set_printoptions(formatter={'float_kind':float_formatter})

# remove the pairs of (r,d) whose d<=0.0
def remove_density_in_profile( r: float, d: float ) -> tuple[float, float]:

    r_no_zero = r[d > 0.0]
    d_no_zero = d[d > 0.0]

    return r_no_zero, d_no_zero

def main() -> None:

    # example density profile, where some values of density is 0.0
    radius  = np.array([  1.0, 2.0,  3.0, 4.0,  5.0, 6.0,  7.0])
    density = np.array([100.0, 0.0, 90.0, 0.0, 80.0, 0.0, 70.0])

    print( f'' )
    print( f'Before:' )
    print( f'{radius                = }' )
    print( f'{density               = }' )

    # remove the zero density points in the density profile
    radius_removed_zero, density_removed_zero = remove_density_in_profile(radius, density)

    print( f'' )
    print( f'After:' )
    print( f'{radius_removed_zero   = }' )
    print( f'{density_removed_zero  = }' )

if __name__ == '__main__':
    main()
