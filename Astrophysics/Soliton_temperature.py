#!/usr/bin/env python3



def Soliton_rho_c_from_r_c(m22: float, r_c: float) -> float:
    """Summary line.

    Extended description of function.

    Parameters
    ----------
    arg1 : int
        Description of arg1
    arg2 : str
        Description of arg2
    *args
        Variable length argument list.
    **kwargs
        Arbitrary keyword arguments.

    Returns
    -------
    bool
        Description of return value

    Raises
    ------
    AttributeError
        The ``Raises`` section is a list of all exceptions
        that are relevant to the interface.
    ValueError
        If `arg1` is equal to `arg2`.

    """
    return 1.945e7/( m22**2 * (r_c*UNIT_L/Const_kpc)**4 )*(Const_Msun/Const_kpc**3)/UNIT_D


def Soliton_fitting_analytical_dens(r, m22, r_c):
    rho_c = Soliton_rho_c_from_r_c(m22, r_c)
    x = r / r_c
    return rho_c*( 1.0+9.06e-2*(x**2) )**(-8)


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

    print('Hello World!')

if __name__ == '__main__':
    main()
