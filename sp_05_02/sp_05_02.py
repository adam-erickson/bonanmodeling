# Supplemental program 5.2

# Author: Adam Erickson, PhD, Washington State University

# ---------------------------------------------------------
# Diurnal cycle of soil temperature with phase change using
# "excess heat" or "apparent heat capacity"
# ---------------------------------------------------------

import numpy as np
from matplotlib import ticker
from matplotlib import pyplot as plt

from soil_thermal_properties import soil_thermal_properties
from soil_temperature import soil_temperature


# --- Physical constant and soil texture variable classes

class PhysCon:
    def __init__(self):
        self.tfrz = 273.15          # Freezing point of water (K)
        self.cwat = 4188.0          # Specific heat of water (J/kg/K)
        self.cice = 2117.27         # Specific heat of ice (J/kg/K)
        self.rhowat = 1000.0        # Density of water (kg/m3)
        self.rhoice = 917.0         # Density of ice (kg/m3)
        self.cvwat = self.cwat * self.rhowat  # Heat capacity of water (J/m3/K)
        self.cvice = self.cice * self.rhoice  # Heat capacity of ice (J/m3/K)
        self.tkwat = 0.57           # Thermal conductivity of water (W/m/K)
        self.tkice = 2.29           # Thermal conductivity of ice (W/m/K)
        self.hfus = 0.3337e6        # Heat of fusion for water at 0 C (J/kg)


class SoilVar:
    """
    Soil texture classes (Cosby et al. 1984. Water Resources Research 20:682-690):
     0: sand
     1: loamy sand
     2: sandy loam
     3: silty loam
     4: loam
     5: sandy clay loam
     6  silty clay loam
     7: clay loam
     8: sandy clay
     9: silty clay
    10: clay
    """

    def __init__(self, soil_texture, method, nsoi, dz):
        # Volumetric soil water content at saturation (porosity)
        # (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)
        self.watsat = np.array(
            [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482])
        # Percent silt
        self.silt = np.array(
            [5.0, 12.0, 32.0, 70.0, 39.0, 15.0, 56.0, 34.0,  6.0, 47.0, 20.0])
        # Percent sand
        self.sand = np.array(
            [92.0, 82.0, 58.0, 17.0, 43.0, 58.0, 10.0, 32.0, 52.0,  6.0, 22.0])
        # Percent clay
        self.clay = np.array(
            [3.0,  6.0, 10.0, 13.0, 18.0, 27.0, 34.0, 34.0, 42.0, 47.0, 58.0])
        # Soil texture class
        self.soil_texture = soil_texture
        # Use excess heat for phase change
        self.method = method
        # Number of layers in soil profile
        self.nsoi = nsoi
        # Soil layer thickness (m)
        self.dz = np.full(shape=nsoi, fill_value=dz)
        # Soil depth (m) at center of layer i (negative distance from surface)
        self.z = np.zeros(shape=nsoi)
        # Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)
        self.z_plus_onehalf = np.zeros(shape=nsoi)
        # Thickness between between z(i) and z(i+1)
        self.dz_plus_onehalf = np.zeros(shape=nsoi)
        # Soil temperature (K)
        self.tsoi = np.zeros(shape=nsoi)
        # Ice fraction of soil water at saturation (kg H2O/m2)
        self.h2osoi_ice = np.zeros(shape=nsoi)
        # Liquid fraction of soil water at saturation (kg H2O/m2)
        self.h2osoi_liq = np.zeros(shape=nsoi)
        # Thermal conducitivty (W/m/K)
        self.tk = np.zeros(shape=nsoi)
        # Volumetric heat capacity (J/m3/K)
        self.cv = np.zeros(shape=nsoi)
        # Energy flux into soil (W/m2)
        self.gsoi = None
        # Soil phase change energy flux (W/m2)
        self.hfsoi = None


# --- Initialize physical constants

physcon = PhysCon()

# --- Initialize soil variables

# sand texture, 'excess-heat' phase-change method, 120 layers in soil profile, soil layer thickness of 0.025 m; 'apparent-heat-capacity' is also available

soilvar = SoilVar(soil_texture=0, method='excess-heat', nsoi=120, dz=0.025)

# --- Model run control parameters

# Mean daily air temperature for diurnal cycle (K)
tmean = physcon.tfrz + 15.0
trange = 10.0                   # Temperature range for diurnal cycle (K)
dt = 1800                       # Time step (seconds)
nday = 200                      # Number of days

# --- Initialize soil layer variables

# Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)

soilvar.z_plus_onehalf[0] = -soilvar.dz[0]

for i in range(1, soilvar.nsoi):
    soilvar.z_plus_onehalf[i] = soilvar.z_plus_onehalf[i-1] - soilvar.dz[i]

# Soil depth (m) at center of layer i (negative distance from surface)

soilvar.z[0] = 0.5 * soilvar.z_plus_onehalf[0]

for i in range(1, soilvar.nsoi):
    soilvar.z[i] = 0.5 * (soilvar.z_plus_onehalf[i-1] +
                          soilvar.z_plus_onehalf[i])

# Thickness between between z(i) and z(i+1)

for i in range(soilvar.nsoi-1):
    soilvar.dz_plus_onehalf[i] = soilvar.z[i] - soilvar.z[i+1]

soilvar.dz_plus_onehalf[soilvar.nsoi-1] = 0.5 * soilvar.dz[soilvar.nsoi-1]

# Initial soil temperature (K) and unfrozen and frozen water (kg H2O/m2)

for i in range(soilvar.nsoi):

    # Temperature (K)

    soilvar.tsoi[i] = physcon.tfrz + 2.0

    # Soil water at saturation (kg H2O/m2)

    h2osoi_sat = soilvar.watsat[soilvar.soil_texture] * \
        physcon.rhowat * soilvar.dz[i]

    # Actual water content is some fraction of saturation

    if soilvar.tsoi[i] > physcon.tfrz:
        soilvar.h2osoi_ice[i] = 0
        soilvar.h2osoi_liq[i] = 0.8 * h2osoi_sat
    else:
        soilvar.h2osoi_ice[i] = 0.8 * h2osoi_sat
        soilvar.h2osoi_liq[i] = 0

# --- Time stepping loop to increment soil temperature

# Counter for output file

m = 0

# Main loop is NTIM iterations per day with a time step of DT seconds.
# This is repeated NDAY times.

ntim = np.int(np.round(86400/dt))

# -- Placeholders

# shape = 1968 = 48 (ntim) * 41 (?)
hour_vec = np.zeros(shape=1968)
z_vec = np.zeros(shape=1968)
tsoi_vec = np.zeros(shape=1968)

hour_out = np.zeros(shape=1968)
z_out = np.zeros(shape=1968)
tsoi_out = np.zeros(shape=1968)

for iday in range(1, nday+1):

    print("day = {:6.0f}".format(iday))

    for itim in range(ntim):

        # Hour of day

        hour = itim * (dt/86400 * 24)

        # Surface temperature: Constant value TMEAN if TRANGE = 0. Otherwise, use a sine
        # wave with max (TMEAN + 1/2 TRANGE) at 2 pm and min (TMEAN - 1/2 TRANGE) at 2 am

        tsurf = tmean + 0.5 * trange * np.sin(2*np.pi/24 * (hour-8.0))

        # Thermal conductivity and heat capacity

        soilvar = soil_thermal_properties(physcon, soilvar)

        # Soil temperatures

        soilvar = soil_temperature(physcon, soilvar, tsurf, dt)

        # Save hourly soil temperature profile for last day

        if iday == nday:

            # Surface output

            # Vector format - to write to data file
            hour_vec[m] = hour
            z_vec[m] = 0
            tsoi_vec[m] = tsurf - physcon.tfrz  # deg C

            # For contour plot
            hour_out[itim] = hour
            z_out[m] = 0
            tsoi_out[m] = tsurf - physcon.tfrz  # deg C

            # Soil layers for top 100 cm

            for i in range(soilvar.nsoi):
                if soilvar.z[i] > -1.0:
                    m += 1
                    hour_vec[m] = hour
                    z_vec[m] = soilvar.z[i] * 100  # cm
                    tsoi_vec[m] = soilvar.tsoi[i] - physcon.tfrz  # deg C
                    z_out[i+1] = soilvar.z[i] * 100  # cm
                    tsoi_out[i+1] = soilvar.tsoi[i] - physcon.tfrz  # deg C

# --- Write output file

a = np.array(["hour", "z", "tsoi"])
header = "{:>12s} {:>12s} {:>12s}".format(*a)

np.savetxt("data_python.txt", np.transpose([hour_vec, z_vec, tsoi_vec]), comments="",
           delimiter=" ", fmt='%12.3f', header=header)

# --- Make contour plot

fig, ax = plt.subplots()
CS = ax.contour(x=hour_out, y=z_out, z=tsoi_out)
ax.clabel(CS, inline=1, fontsize=10)
ax.set_title("Simplest default with labels")
plt.show()

# Replicate Matlab style
# ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
# ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax.set_xlim(0, 1)
#ax.set_ylim(0, 3)
#ax.tick_params(axis='both', direction='in', top=True, right=True)
# ax.legend(loc='best', fontsize='small', edgecolor='k',
#          fancybox=False, framealpha=1, borderaxespad=1)
#plt.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.93)

#plt.title(r"Soil Temperature ($\mathdefault{^oC}$)")
#plt.xlabel("Time of day (hours)")
#plt.ylabel("Soil depth (cm)")
# plt.show()
