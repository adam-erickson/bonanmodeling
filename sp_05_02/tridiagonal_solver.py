# Author: Adam Erickson, PhD, Washington State University

import numpy as np


def tridiagonal_solver(a, b, c, d, n):
    """
    Solve for U given the set of equations R * U = D, where U is a vector
    of length N, D is a vector of length N, and R is an N x N tridiagonal
    matrix defined by the vectors A, B, C each of length N. A(1) and
    C(N) are undefined and are not referenced.

            |B(1) C(1) ...  ...  ...                     |
            |A(2) B(2) C(2) ...  ...                     |
        R = |     A(3) B(3) C(3) ...                     |
            |                    ... A(N-1) B(N-1) C(N-1)|
            |                    ... ...    A(N)   B(N)  |

    The system of equations is written as:

        A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i

    for i = 1 to N. The solution is found by rewriting the equations
    so that:

        U_i = F_i - E_i * U_i+1
    """

    e = np.zeros(n)
    f = np.zeros(n)
    u = np.zeros(n)

    # --- Forward sweep (1 -> N) to get E and F

    e[0] = c[0] / b[0]

    for i in range(1, n-1, 1):
        e[i] = c[i] / (b[i] - a[i] * e[i-1])

    f[0] = d[0] / b[0]

    for i in range(1, n, 1):
        f[i] = (d[i] - a[i] * f[i-1]) / (b[i] - a[i] * e[i-1])

    # --- Backward substitution (N -> 1) to solve for U

    u[n] = f[n]

    for i in range(n-1, 1, -1):
        u[i] = f[i] - e[i] * u[i+1]

    return u
