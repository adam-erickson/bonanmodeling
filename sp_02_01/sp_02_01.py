# Supplemental program 2.1

# Author: Adam Erickson, PhD, Washington State University

# ----------------------------------------------
# Calculate and graph leaf area density profiles
# ----------------------------------------------

import numpy as np
from scipy.special import beta
from matplotlib import pyplot as plt

# --- Parameters for beta distribution

p = np.array([2.5, 3.5, 11.5])
q = np.array([2.5, 2.0,  3.5])

# --- Canopy parameters

LAI = 5              # leaf area index (m2/m2)
hc = 10              # canopy height (m)

print("Leaf area index = {:6.2f}".format(LAI))

# --- Create a vector of heights (z) with linearly spaced
# values from z_min to z_max with increment dz

z_min = 0            # minimum
z_max = hc           # maximum
dz = 0.1             # increment
# z is linearly spaced between z_min and z_max with increment dz
z = np.arange(start=z_min, stop=z_max+dz, step=dz, dtype=np.float64)

# --- Parameters for plotting

y1 = np.zeros(len(z))
y2 = np.zeros(len(z))
y3 = np.zeros(len(z))

# --- Calculate the leaf area density profile for each [p,q]

# Loop over each p,q
for i in range(3):

    # Loop over each height
    sigma = 0
    for j in range(len(z)):
        x = z[j] / hc
        lad = (LAI / hc) * (x ** (p[i]-1) * (1 - x)
                            ** (q[i]-1)) / beta(p[i], q[i])

        # Numerically sum leaf area for each height
        sigma += lad * dz

        # Save output for graphing
        if i == 0:
            y1[j] = lad
        elif i == 1:
            y2[j] = lad
        elif i == 2:
            y3[j] = lad
        else:
            print("Index outside range of z")

    print("p, q = {:6.2f} {:6.2f}".format(p[i], q[i]))
    print("Leaf area index (numerical) = {:8.4f}".format(sigma))

# --- Make graph for leaf area density in relation to relative height (z/hc)

z_rel = z / hc

fig, ax = plt.subplots()
ax.plot(y1, z_rel, 'b-', label="p,q = 2.5,2.5")
ax.plot(y2, z_rel, 'r-', label="p,q = 3.5,2.0")
ax.plot(y3, z_rel, 'g-', label="p,q = 11.5,3.5")

ax.tick_params(axis='both', direction='in', top=True, right=True)
ax.legend(loc='lower right', fontsize='small', framealpha=1)

plt.title("Profiles", fontweight='bold')
plt.xlabel(r"Leaf area density ($\mathdefault{m^2 m^{-3}}$)")
plt.ylabel(r"Height ($\mathdefault{z/h_c}$)")
plt.show()
