#!/usr/bin/env python3

import numpy as np

def load_from_table( Filename_table: str ) -> tuple[float, float]:
    """This funcion load a two-column text file and return the two colums as two arrays.

    Parameters
    ----------
    Filename_table : str
        The filename of the input text file

    Returns
    -------
    Table_x : float
        The array for the fisrt column
    Table_y : float
        The array for the second colum

    """

    # load the table
    Table_x, Table_y = np.loadtxt( Filename_table, skiprows=1, unpack=True )

    return Table_x, Table_y


def save_table_to_file( Filename_output: str, Table_x: float, Table_y: float ) -> None:
    """This function save two arrays and save them into two columns in a text file.

    Parameters:
    ----------
    Filename_output : str
        The filename of the output text file
    Table_x         : float
        The array for the fisrt column
    Table_y         : float
        Yhe array for the second colum

    Returns
    -------
    None

    """

    # save the table to text file
    np.savetxt( Filename_output,
                np.column_stack( (Table_x, Table_y) ),
                fmt = '          %9.8e',
                header = '         x (unit_of_x)            y (unit_of_y)' )

    return

def main() -> None:
    """This function is a example to read two columns from the text file and save them into another text file.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    # readi a density profile
    radius, density = load_from_table( "Density_table_input" )

    # save the density profile to file
    save_table_to_file( "Density_table_output", radius, density )

    return

if __name__ == '__main__':
    main()
