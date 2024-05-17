#!/usr/bin/env python3

import numpy as np

def main() -> None:
    """main function
    This function is a template.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    NEWTON_G  = 1.0000000e+00
    Total_M   = 2.3062183e-02
    Max_R     = 1.5000000e+00

    Ave_Rho   = 3.0*Total_M/(4.0*np.pi*Max_R**3)
    t_ff      = ( 3.0*np.pi/(32.0*NEWTON_G*Ave_Rho) )**0.5

    print( f'{NEWTON_G =: .8e}' )
    print( f'{Total_M  =: .8e}' )
    print( f'{Max_R    =: .8e}' )
    print( f'{Ave_Rho  =: .8e}' )
    print( f'{t_ff     =: .8e}' )

if __name__ == '__main__':
    main()


