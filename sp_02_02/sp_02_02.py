# Supplemental program 2.2

# ------------------------------------------------------------------------------
# Evaluate the leaf angle probability density function (PDF) from the beta
# distribution using the mean and standard deviation of the leaf inclination
# angle. The leaf angle PDF is calculated for 9 angle classes between 5 and 85
# degrees in increments of 10 degrees.
# ------------------------------------------------------------------------------

import numpy as np
from scipy.special import beta
import scipy.integrate as integrate
from matplotlib import ticker
from matplotlib import pyplot as plt

# --- The variable "leaf" defines the leaf angle distribution type with parameters:
# lad_ave - mean leaf angle (radians)
# lad_std - standard deviation of leaf angle (radians)

leaf = 'Planophile'
# leaf = 'Erectophile';
# leaf = 'Plagiophile';
# leaf = 'Uniform';
# leaf = 'Spherical';

if leaf == 'Planophile':
    lad_ave = 26.76 * (np.pi/180)
    lad_std = 18.5068 * (np.pi/180)
elif leaf == 'Erectophile':
    lad_ave = 63.24 * (np.pi/180)
    lad_std = 18.4960 * (np.pi/180)
elif leaf == 'Plagiophile':
    lad_ave = 45.00 * (np.pi/180)
    lad_std = 16.2681 * (np.pi/180)
elif leaf == 'Uniform':
    lad_ave = 45.00 * (np.pi/180)
    lad_std = 25.9808 * (np.pi/180)
elif leaf == 'Spherical':
    lad_ave = 57.30 * (np.pi/180)
    lad_std = 21.5485 * (np.pi/180)
else:
    "Invalid leaf type"

# --- Convert these to the p,q parameters for the beta distribution

num = 1 - (lad_std*lad_std + lad_ave*lad_ave) / (lad_ave * np.pi / 2)
den = (lad_std*lad_std + lad_ave*lad_ave) / (lad_ave*lad_ave) - 1
p = num / den
q = ((np.pi/2) / lad_ave - 1) * p

# --- Calculate leaf inclination angle probability density function (PDF) and
# fractional abundance (lad) for 9 10-degree bins

# Leaf inclination angle increment (radians)
dangle = 10 * (np.pi/180)
# Leaf inclination angle (degrees)
angle = np.array([5, 15, 25, 35, 45, 55, 65, 75, 85])
angle = np.multiply(angle, (np.pi/180))    # degrees -> radians

# --- Loop through each angle
beta_pdf = np.zeros(len(angle))
beta_lad = np.zeros(len(angle))

for i in range(len(angle)):
    x = angle[i] / (np.pi/2)
    fp = x ** (p - 1)
    fq = (1 - x) ** (q - 1)
    # Leaf angle probability density function
    beta_pdf[i] = 2 / np.pi * fp * fq / beta(p, q)
    # Fraction of leaves in this angle bin
    beta_lad[i] = beta_pdf[i] * dangle

# --- Calculate the known solution
exact_pdf = np.zeros(len(angle))
exact_lad = np.zeros(len(angle))

for i in range(len(angle)):

    # Exact leaf angle probability density function
    if leaf == 'Planophile':
        exact_pdf[i] = 2 / np.pi * (1 + np.cos(2 * angle[i]))
    elif leaf == 'Erectophile':
        exact_pdf[i] = 2 / np.pi * (1 - np.cos(2 * angle[i]))
    elif leaf == 'Plagiophile':
        exact_pdf[i] = 2 / np.pi * (1 - np.cos(4 * angle[i]))
    elif leaf == 'Uniform':
        exact_pdf[i] = 2 / np.pi
    elif leaf == 'Spherical':
        exact_pdf[i] = np.sin(angle[i])
    else:
        "Invalid leaf type"

    # Exact relative leaf angle distribution (fraction)
    exact_lad[i] = exact_pdf[i] * dangle

# --- Print out fractional abundance and compare with known solution
# sum -  sum of PDF * dangle (this sums to 1)
# ave -  sum of angle * PDF * dangle (this is the mean)

beta_sum = 0
beta_ave = 0
exact_sum = 0
exact_ave = 0

print("Leaf type = {:15s}".format(leaf))
print("     Angle         beta            exact ")

for i in range(len(angle)):
    beta_sum = beta_sum + beta_lad[i]
    beta_ave = beta_ave + angle[i] * beta_lad[i]
    exact_sum = exact_sum + exact_lad[i]
    exact_ave = exact_ave + angle[i] * exact_lad[i]
    print("{:10.2f} {:15.4f} {:15.4f}".format(angle[
        i]*180/np.pi, beta_lad[i], exact_lad[i]))

print("beta distribution")
print("Sum of leaf angle distribution = {:15.4f}".format(beta_sum))
print("Mean leaf angle = {:15.4f}".format(beta_ave*180/np.pi))
print("Exact solution")
print("Sum of leaf angle distribution = {:15.4f}".format(exact_sum))
print("Mean leaf angle = {:15.4f}".format(exact_ave*180/np.pi))

# --- Analytical mean leaf angle

if leaf == 'Planophile':
    def fx(x): return x * 2 / np.pi * (1 + np.cos(2 * x))
elif leaf == 'Erectophile':
    def fx(x): return x * 2 / np.pi * (1 - np.cos(2 * x))
elif leaf == 'Plagiophile':
    def fx(x): return x * 2 / np.pi * (1 - np.cos(4 * x))
elif leaf == 'Uniform':
    def fx(x): return x * 2 / np.pi
elif leaf == 'Spherical':
    def fx(x): return x * np.sin(x)
else:
    "Invalid leaf type"

analytical_ave, err = integrate.quad(fx, 0, np.pi/2)

print("Analytical solution")
print("Mean leaf angle = {:15.4f}".format(analytical_ave*180/np.pi))

# --- Calculate Ross index

F1 = beta_lad[0] + beta_lad[1] + beta_lad[2]
F2 = beta_lad[3] + beta_lad[4] + beta_lad[5]
F3 = beta_lad[6] + beta_lad[7] + beta_lad[8]
beta_xl = 0.5 * (np.abs(0.134-F1) + np.abs(0.366-F2) + np.abs(0.5-F3))
if (0.5-F3) < 0:
    beta_xl = -beta_xl

F1 = exact_lad[0] + exact_lad[1] + exact_lad[2]
F2 = exact_lad[3] + exact_lad[4] + exact_lad[5]
F3 = exact_lad[6] + exact_lad[7] + exact_lad[8]
exact_xl = 0.5 * (np.abs(0.134-F1) + np.abs(0.366-F2) + np.abs(0.5-F3))
if (0.5-F3) < 0:
    exact_xl = -exact_xl

print("Ross index")
print("beta distribution = {:15.4f}".format(beta_xl))
print("   Exact solution = {:15.4f}".format(exact_xl))

# --- Graph PDFs

# x is linearly spaced values between 0 and pi/2 with increment (pi/2)/100
dx = (np.pi/2) / 100
x = np.arange(start=0, stop=(np.pi/2)+dx, step=dx)

if leaf == 'Planophile':
    y = 2 / np.pi * (1 + np.cos(2*x))
elif leaf == 'Erectophile':
    y = 2 / np.pi * (1 - np.cos(2*x))
elif leaf == 'Plagiophile':
    y = 2 / np.pi * (1 - np.cos(4*x))
elif leaf == 'Uniform':
    y = 2 / np.pi
elif leaf == 'Spherical':
    y = np.sin(x)
else:
    "Invalid leaf type"

x = np.multiply(x, (180/np.pi))           # radians -> degrees
angle = np.multiply(angle, (180/np.pi))   # radians -> degrees

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(angle, beta_pdf, 'b--o', linewidth=1, label="beta")
ax.plot(angle, exact_pdf, 'r--*', linewidth=1, label="exact")
ax.plot(x, y, 'g', linewidth=1, label="analytical")

# Replicate Matlab style
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.set_xlim(0, 90)
ax.set_ylim(0, 1.4)
ax.tick_params(axis='both', direction='in', top=True, right=True)
ax.legend(loc='best', fontsize='small', edgecolor='k',
          fancybox=False, framealpha=1, borderaxespad=1)
plt.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.93)

plt.title(leaf, fontweight='bold', fontdict={'fontsize': 10})
plt.xlabel("Leaf angle (degrees)")
plt.ylabel("PDF")
plt.show()
