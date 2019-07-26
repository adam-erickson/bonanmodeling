# Supplemental program 5.1

# Author: Adam Erickson, PhD, Washington State University

# ----------------------------------------
# Calculate and graph thermal conductivity
# ----------------------------------------

import numpy as np
from matplotlib import ticker
from matplotlib import pyplot as plt

# Constants

cwat = 4188               # Specific heat of water (J/kg/K)
cice = 2117.27            # Specific heat ice (J/kg/K)

rho_wat = 1000            # Density of water (kg/m3)
rho_ice = 917             # Density of ice (kg/m3)

cvwat = cwat * rho_wat    # Heat capacity of water (J/m3/K)
cvice = cice * rho_ice    # Heat capacity of ice (J/m3/K)
cvsol = 1.926e06          # Heat capacity of soil solids (J/m3/K)

tkwat = 0.57              # Thermal conductivity of water (W/m/K)
tkice = 2.29              # Thermal conductivity of ice (W/m/K)
tk_quartz = 7.7           # Thermal conductivity of quartz (W/m/K)

# Soil texture classes (Cosby et al. 1984. Water Resources Research 20:682-690)

#  1: sand
#  2: loamy sand
#  3: sandy loam
#  4: silty loam
#  5: loam
#  6: sandy clay loam
#  7  silty clay loam
#  8: clay loam
#  9: sandy clay
# 10: silty clay
# 11: clay

silt = np.array([5, 12, 32, 70, 39, 15, 56, 34,  6, 47, 20])    # Percent silt
sand = np.array([92, 82, 58, 17, 43, 58, 10, 32, 52,  6, 22])   # Percent sand
clay = np.array([3,  6, 10, 13, 18, 27, 34, 34, 42, 47, 58])    # Percent clay

# Volumetric soil water content at saturation (porosity)
# (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)

watsat = np.array([0.395, 0.410, 0.435, 0.485, 0.451,
                   0.420, 0.477, 0.476, 0.426, 0.492, 0.482])

# Define 5 soil types to process

soiltyp = np.array([0, 2, 4, 7, 10])

# Set relative soil water content (s) from 0 to 1

inc = 0.05                       # increment
n = np.int((1 - 0) / inc + 1)    # number of values
s = np.linspace(0, 1, n)         # n evenly spaced values between 0 and 1

# Placeholder arrays for results

tk1 = np.zeros(len(s))
cv1 = np.zeros(len(s))
tk2 = np.zeros(len(s))
cv2 = np.zeros(len(s))
tk3 = np.zeros(len(s))
cv3 = np.zeros(len(s))
tk4 = np.zeros(len(s))
cv4 = np.zeros(len(s))
tk5 = np.zeros(len(s))
cv5 = np.zeros(len(s))

# Loop through each soil type

for i in range(len(soiltyp)):

    # Soil texture to process

    k = soiltyp[i]

    # Thermal conductivity and heat capacity for each soil moisture

    for j in range(len(s)):

        # Volumetric water content

        h2osoi = s[j] * watsat[k]

        # Dry thermal conductivity (W/m/K) from bulk density (kg/m3)

        bd = 2700 * (1 - watsat[k])
        tkdry = (0.135 * bd + 64.7) / (2700 - 0.947 * bd)

        # Soil solids thermal conducitivty (W/m/K) from quartz fraction
        # tko = thermal conductivity of other minerals (W/m/K)

        quartz = sand[k] / 100
        if quartz > 0.2:
            tko = 2.0
        else:
            tko = 3.0

        tksol = (tk_quartz ** quartz) * (tko ** (1-quartz))

        # Unfrozen and frozen saturated thermal conductivity (W/m/K)

        tksat_u = (tksol**(1-watsat[k])) * (tkwat**watsat[k])
        tksat_f = (tksol**(1-watsat[k])) * (tkice**watsat[k])

        # Unfrozen and frozen Kersten number

        if sand[k] < 50:
            ke_u = np.log10(np.max([s[j], 0.1])) + 1
        else:
            ke_u = 0.7 * np.log10(np.max([s[j], 0.05])) + 1

        ke_f = s[j]

        # Unfrozen and frozen thermal conductivity (W/m/K)

        tku = (tksat_u - tkdry) * ke_u + tkdry
        tkf = (tksat_f - tkdry) * ke_f + tkdry

        # Unfrozen and frozen heat capacity (J/m3/K)

        cvu = (1 - watsat[k]) * cvsol + cvwat * h2osoi
        cvf = (1 - watsat[k]) * cvsol + cvice * h2osoi

        # Save values for each texture type

        if i == 0:
            tk1[j] = tku
            cv1[j] = cvu * 1e-06
        elif i == 1:
            tk2[j] = tku
            cv2[j] = cvu * 1e-06
        elif i == 2:
            tk3[j] = tku
            cv3[j] = cvu * 1e-06
        elif i == 3:
            tk4[j] = tku
            cv4[j] = cvu * 1e-06
        elif i == 4:
            tk5[j] = tku
            cv5[j] = cvu * 1e-06
        else:
            print("Error: value i out of range")

# Make graph

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(s, tk1, 'r-', linewidth=1, label="sand")
ax.plot(s, tk2, 'g-', linewidth=1, label="sandy loam")
ax.plot(s, tk3, 'b-', linewidth=1, label="loam")
ax.plot(s, tk4, 'm-', linewidth=1, label="clay loam")
ax.plot(s, tk5, 'c-', linewidth=1, label="clay")

# Replicate Matlab style
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.set_xlim(0, 1)
ax.set_ylim(0, 3)
ax.tick_params(axis='both', direction='in', top=True, right=True)
ax.legend(loc='best', fontsize='small', edgecolor='k',
          fancybox=False, framealpha=1, borderaxespad=1)
plt.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.93)

plt.title("Thermal conductivity", fontweight='bold', fontdict={'fontsize': 10})
plt.xlabel("Relative soil moisture")
plt.ylabel(r"Thermal conductivity ($\mathdefault{W m^{-1} K^{-1}}$)")
plt.show()

# Write formated output to text file: n rows x 6 columns
# column 1 = relative soil water (s)
# column 2 = thermal conductivity for soil 1 (W/m/K)
# column 3 = thermal conductivity for soil 2 (W/m/K)
# column 4 = thermal conductivity for soil 3 (W/m/K)
# column 5 = thermal conductivity for soil 4 (W/m/K)
# column 6 = thermal conductivity for soil 5 (W/m/K)

# First, create a header string with proper spacing
a = np.array(["s", "tex1", "tex2", "tex3", "tex4", "tex5"])
header = "{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}".format(*a)

np.savetxt("tk_python.txt", np.transpose([s, tk1, tk2, tk3, tk4, tk5]), comments="",
           delimiter=" ", fmt='%12.4f', header=header)

# Or, you can write to the file manually
# with open("tk_python.txt", 'w') as f:
#    for i in ["s", "tex1", "tex2", "tex3", "tex4", "tex5"]:
#        f.write("{:>12s}".format(i))
#    f.write("\n")
#    for j in np.transpose([s, tk1, tk2, tk3, tk4, tk5]):
#        for k in j:
#            f.write("{:>12.4f}".format(k))
#        f.write("\n")

# Write formated output to text file: n rows x 6 columns
# column 1 = relative soil water (s)
# column 2 = heat capacity for soil 1 (MJ/m3/K)
# column 3 = heat capacity for soil 2 (MJ/m3/K)
# column 4 = heat capacity for soil 3 (MJ/m3/K)
# column 5 = heat capacity for soil 4 (MJ/m3/K)
# column 6 = heat capacity for soil 5 (MJ/m3/K)

np.savetxt("cv_python.txt", np.transpose([s, cv1, cv2, cv3, cv4, cv5]), comments="",
           delimiter=" ", fmt='%12.4f', header=header)
