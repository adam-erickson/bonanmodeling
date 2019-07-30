# Author: Adam Erickson, PhD, Washington State University

import numpy as np

from tridiagonal_solver import tridiagonal_solver
from phase_change import phase_change


def soil_temperature(physcon, soilvar, tsurf, dt):
    """
    Use an implicit formulation with the surface boundary condition specified as the surface temperature to solve for soil temperatures at time n+1.

    Calculate soil temperatures as:

          dT   d     dT
       cv -- = -- (k --)
          dt   dz    dz

    where: T = temperature (K)
            t = time (s)
            z = depth (m)
            cv = volumetric heat capacity (J/m3/K)
            k = thermal conductivity (W/m/K)

    Set up a tridiagonal system of equations to solve for T at time n+1, where the temperature equation for layer i is

       d_i = a_i [T_i-1] n+1 + b_i [T_i] n+1 + c_i [T_i+1] n+1

    For soil layers undergoing phase change, set T_i = Tf (freezing) and use excess energy to freeze or melt ice:

       Hf_i = (Tf - [T_i] n+1) * cv_i * dz_i / dt

    During the phase change, the unfrozen and frozen soil water (h2osoi_liq, h2osoi_ice) are adjusted.

    Or alternatively, use the apparent heat capacity method to account for phase change. In this approach, h2osoi_liq and h2osoi_ice are not calculated.

    ------------------------------------------------------
    # Input
       tsurf                   ! Surface temperature (K)
       dt                      ! Time step (s)
       soilvar.method          ! Use excess heat or apparent heat capacity for phase change
       soilvar.nsoi            ! Number of soil layers
       soilvar.z               ! Soil depth (m)
       soilvar.z_plus_onehalf  ! Soil depth (m) at i+1/2 interface between layers i and i+1
       soilvar.dz              ! Soil layer thickness (m)
       soilvar.dz_plus_onehalf ! Thickness (m) between between i and i+1
       soilvar.tk              ! Thermal conductivity (W/m/K)
       soilvar.cv              ! Heat capacity (J/m3/K)

    # Input/output
       soilvar.tsoi            ! Soil temperature (K)
       soilvar.h2osoi_liq      ! Unfrozen water, liquid (kg H2O/m2)
       soilvar.h2osoi_ice      ! Frozen water, ice (kg H2O/m2)

    # Output
       soilvar.gsoi            ! Energy flux into soil (W/m2)
       soilvar.hfsoi           ! Soil phase change energy flux (W/m2)
    ------------------------------------------------------
    """

    """
    Matlab tsurf
    283.530602

    Python tsurf
    283.819873
    """

    # --- Placeholder arrays

    tk_plus_onehalf = np.zeros(soilvar.nsoi)
    a = np.zeros(soilvar.nsoi)
    c = np.zeros(soilvar.nsoi)
    b = np.zeros(soilvar.nsoi)
    d = np.zeros(soilvar.nsoi)

    # --- Save current soil temperature for energy conservation check

    tsoi0 = soilvar.tsoi

    # --- Thermal conductivity at interface (W/m/K)

    for i in range(soilvar.nsoi-1):
        tk_plus_onehalf[i] = soilvar.tk[i] * soilvar.tk[i+1] * (soilvar.z[i]-soilvar.z[i+1]) / \
            (soilvar.tk[i]*(soilvar.z_plus_onehalf[i]-soilvar.z[i+1]) +
             soilvar.tk[i+1]*(soilvar.z[i]-soilvar.z_plus_onehalf[i]))

    # --- Set up tridiagonal matrix

    # Top soil layer with tsurf as boundary condition

    i = 0
    m = soilvar.cv[i] * soilvar.dz[i] / dt
    a[i] = 0
    c[i] = -tk_plus_onehalf[i] / soilvar.dz_plus_onehalf[i]
    b[i] = m - c[i] + soilvar.tk[i] / (0 - soilvar.z[i])
    d[i] = m * soilvar.tsoi[i] + soilvar.tk[i] / (0 - soilvar.z[i]) * tsurf

    # Layers 2 to nsoi-1

    for i in range(1, soilvar.nsoi-1):
        m = soilvar.cv[i] * soilvar.dz[i] / dt
        a[i] = -tk_plus_onehalf[i-1] / soilvar.dz_plus_onehalf[i-1]
        c[i] = -tk_plus_onehalf[i] / soilvar.dz_plus_onehalf[i]
        b[i] = m - a[i] - c[i]
        d[i] = m * soilvar.tsoi[i]

    # Bottom soil layer with zero heat flux

    i = soilvar.nsoi-1
    m = soilvar.cv[i] * soilvar.dz[i] / dt
    a[i] = -tk_plus_onehalf[i-1] / soilvar.dz_plus_onehalf[i-1]
    c[i] = 0
    b[i] = m - a[i]
    d[i] = m * soilvar.tsoi[i]

    # --- Solve for soil temperature

    soilvar.tsoi = tridiagonal_solver(a, b, c, d, soilvar.nsoi)

    # --- Derive energy flux into soil (W/m2)

    soilvar.gsoi = soilvar.tk[0] * \
        (tsurf - soilvar.tsoi[0]) / (0 - soilvar.z[0])

    # --- Phase change for soil layers undergoing freezing of thawing

    if soilvar.method == 'apparent-heat-capacity':

        # No explicit phase change energy flux. This is included in the heat capacity.

        soilvar.hfsoi = 0

    elif soilvar.method == 'excess-heat':

        # Adjust temperatures for phase change. Freeze or melt ice using energy
        # excess or deficit needed to change temperature to the freezing point.
        # The variable hfsoi is returned as the energy flux from phase change (W/m2).

        soilvar = phase_change(physcon, soilvar, dt)

    # --- Check for energy conservation

    # Sum change in energy (W/m2)

    edif = 0
    for i in range(soilvar.nsoi):
        edif = edif + soilvar.cv[i] * soilvar.dz[i] * \
            (soilvar.tsoi[i] - tsoi0[i]) / dt

    # Error check

    err = edif - soilvar.gsoi - soilvar.hfsoi
    #print("Original values: ", 464.621857, 464.621857, 0.000000)
    #print("Error check values: ", edif, soilvar.gsoi, soilvar.hfsoi)
    if np.abs(err) > 1e-03:
        raise Exception("Soil temperature energy conservation error")

    return soilvar
