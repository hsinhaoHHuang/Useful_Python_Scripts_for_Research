#!/usr/bin/env python3

import numpy as np
import random
from matplotlib import pyplot as plt

def LinearRegression(x: list[float], y: list[float]) -> list[float]:
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
    N = len(x)
    [SX2, SX, SXY, SY] = [sum(x**2), sum(x), sum(x*y), sum(y)]
    [M, V] = [ np.array([ [SX2, SX], [SX, N]]), np.array([ [SXY], [SY] ]) ]
    [a, b] = np.dot( np.linalg.inv(M), V )

    return [a, b]

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
    [a_, b_, N, R] = [2, 10, 100, 6]
    x = np.linspace( 0, 10, N )
    y = a_*x + b_
    for i in range(N):
       y[i] = y[i] + R*(random.random()-0.5)

    [a, b] = LinearRegression( x, y )
    y_f = a*x+b

    fig = plt.figure(figsize=(8,4.5))
    ax = fig.add_subplot(111)
    ax.scatter(x,y)
    ax.plot(x, y_f)
    for i in range(N):
       ax.plot( [x[i], x[i]], [y[i], y_f[i]], ':', color='r', zorder=0 )

    Infostr = 'Fit: y = '+'{:.5f}'.format(float(a))+'x + '+'{:.5f}'.format(float(b))
    ax.set_xlabel('x\n'+Infostr)
    ax.set_ylabel('y')
    ax.set_title('Data points & Line fit')
    ax.legend(['${x_i,y_i}$ Data', 'Line fit', 'Residuals'])

    plt.tight_layout()
    fig.savefig( 'fig_linear_regression.png' )

if __name__ == '__main__':
    main()
