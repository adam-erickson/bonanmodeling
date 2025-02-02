# Author: Adam Erickson, PhD, Washington State University

import numpy as np


def phase_change(physcon, soilvar, dt):
    """
    Adjust temperatures for phase change. Freeze or melt ice using
    energy excess or deficit needed to change temperature to the
    freezing point.

    ------------------------------------------------------
    # Input
      dt                      ! Time step (s)
      physcon.hfus            ! Heat of fusion for water at 0 C (J/kg)
      physcon.tfrz            ! Freezing point of water (K)
      soilvar.nsoi            ! Number of soil layers
      soilvar.dz              ! Soil layer thickness (m)
      soilvar.cv              ! Volumetric heat capacity (J/m3/K)

    # Input/output
      soilvar.tsoi            ! Soil temperature (K)
      soilvar.h2osoi_liq      ! Unfrozen water, liquid (kg H2O/m2)
      soilvar.h2osoi_ice      ! Frozen water, ice (kg H2O/m2)

    # Output
      soilvar.hfsoi           ! Soil phase change energy flux (W/m2)
    ------------------------------------------------------
    """

    # --- Initialize total soil heat of fusion to zero

    soilvar.hfsoi = 0

    # --- Now loop over all soil layers to calculate phase change

    for i in range(soilvar.nsoi):

        # --- Save variables prior to phase change

        # Amount of liquid water before phase change
        wliq0 = soilvar.h2osoi_liq[i]
        wice0 = soilvar.h2osoi_ice[i]     # Amount of ice before phase change
        wmass0 = wliq0 + wice0            # Amount of total water before phase change
        # Soil temperature before phase change
        tsoi0 = soilvar.tsoi[i]

        # --- Identify melting or freezing layers and set temperature to freezing

        # Default condition is no phase change (imelt = 0)

        imelt = 0

        # Melting: if ice exists above melt point, melt some to liquid.
        # Identify melting by imelt = 1

        if soilvar.h2osoi_ice[i] > 0 and soilvar.tsoi[i] > physcon.tfrz:
            imelt = 1
            soilvar.tsoi[i] = physcon.tfrz

        # Freezing: if liquid exists below melt point, freeze some to ice.
        # Identify freezing by imelt = 2

        if soilvar.h2osoi_liq[i] > 0 and soilvar.tsoi[i] < physcon.tfrz:
            imelt = 2
            soilvar.tsoi[i] = physcon.tfrz

        # --- Calculate energy for freezing or melting

        # The energy for freezing or melting (W/m2) is assessed from the energy
        # excess or deficit needed to change temperature to the freezing point.
        # This is a potential energy flux, because cannot melt more ice than is
        # present or freeze more liquid water than is present.
        #
        # heat_flux_pot > 0: freezing; heat_flux_pot < 0: melting

        if imelt > 0:
            heat_flux_pot = (soilvar.tsoi[i] - tsoi0) * \
                soilvar.cv[i] * soilvar.dz[i] / dt
        else:
            heat_flux_pot = 0

        # Maximum energy for melting or freezing (W/m2)

        if imelt == 1:
            heat_flux_max = -soilvar.h2osoi_ice[i] * physcon.hfus / dt
        elif imelt == 2:
            heat_flux_max = soilvar.h2osoi_liq[i] * physcon.hfus / dt

        # --- Now freeze or melt ice

        if imelt > 0:

            # Change in ice (kg H2O/m2/s): freeze (+) or melt (-)

            ice_flux = heat_flux_pot / physcon.hfus

            # Update ice (kg H2O/m2)

            soilvar.h2osoi_ice[i] = wice0 + ice_flux * dt

            # Cannot melt more ice than is present

            soilvar.h2osoi_ice[i] = np.max(
                np.array([0, soilvar.h2osoi_ice[i]]))

            # Ice cannot exceed total water that is present

            soilvar.h2osoi_ice[i] = np.min(
                np.array([wmass0, soilvar.h2osoi_ice[i]]))

            # Update liquid water (kg H2O/m2) for change in ice

            soilvar.h2osoi_liq[i] = np.max(
                np.array([0, wmass0-soilvar.h2osoi_ice[i]]))

            # Actual energy flux from phase change (W/m2). This is equal to
            # heat_flux_pot except if tried to melt too much ice.

            heat_flux = physcon.hfus * (soilvar.h2osoi_ice[i] - wice0) / dt

            # Sum energy flux from phase change (W/m2)

            soilvar.hfsoi = soilvar.hfsoi + heat_flux

            # Residual energy not used in phase change is added to soil temperature

            residual = heat_flux_pot - heat_flux
            soilvar.tsoi[i] = soilvar.tsoi[i] - residual * \
                dt / (soilvar.cv[i] * soilvar.dz[i])

            # Error check: make sure actual phase change does not exceed permissible phase change

            if np.abs(heat_flux) > np.abs(heat_flux_max):
                raise Exception(
                    "Soil temperature energy conservation error: phase change")

            # Freezing: make sure actual phase change does not exceed permissible phase change
            # and that the change in ice does not exceed permissible change

            if imelt == 2:

                # Energy flux (W/m2)

                constraint = np.min(np.array([heat_flux_pot, heat_flux_max]))
                err = heat_flux - constraint
                if np.abs(err) > 1e-03:
                    raise Exception(
                        "Soil temperature energy conservation error: freezing energy flux")

                # Change in ice (kg H2O/m2)

                err = (soilvar.h2osoi_ice(i) - wice0) - \
                    constraint / physcon.hfus * dt
                if np.abs(err) > 1e-03:
                    raise Exception(
                        "Soil temperature energy conservation error: freezing ice flux")

            # Thawing: make sure actual phase change does not exceed permissible phase change
            # and that the change in ice does not exceed permissible change

            elif imelt == 1:

                # Energy flux (W/m2)

                constraint = np.max(np.array([heat_flux_pot, heat_flux_max]))
                err = heat_flux - constraint
                if np.abs(err) > 1e-03:
                    raise Exception(
                        "Soil temperature energy conservation error: thawing energy flux")

                # Change in ice (kg H2O/m2)

                err = (soilvar.h2osoi_ice[i] - wice0) - \
                    constraint / physcon.hfus * dt
                if np.abs(err) > 1e-03:
                    raise Exception(
                        "Soil temperature energy conservation error: thawing ice flux")
    return soilvar
