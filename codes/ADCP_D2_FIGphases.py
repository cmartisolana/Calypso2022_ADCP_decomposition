# Phase diagram of Helmholtz-decomposed ADCP spectra
# Author: Cristina Martí-Solana
# Description:
#   This script visualizes the relationship between rotational (Kψ) and divergent (Kφ) kinetic energy spectra
#   for different layers and eddy phases, using results from the ADCP spectral analysis.
#   Each subplot corresponds to a different layer; colors represent wavelength (1/wavenumber).
#   Intended for publication and reproducibility.

import xarray as xr
import matplotlib.pyplot as plt
import cmocean as cm
import scienceplots
plt.style.use(['science','no-latex',"ieee"])
plt.rcParams['figure.constrained_layout.use'] = True

DIRECTIONS = ["west","east"]

fig, ax = plt.subplots(1, 3, figsize=(4.2, 2), sharex=True, sharey=True)

for dir in DIRECTIONS:
    D1path = f"data/ADCP/D1_files/{dir}.nc"
    D1file = xr.open_dataset(D1path)

    for i in range(3):  # Loop over layers
        # Scatter plot for each phase, colored by wavelength (1/wavenumber)
        ax[i].scatter(D1file.Kpsi[i, 0], D1file.Kphi[i, 0], c=1/D1file.wavenumber, s=6,
                      cmap=cm.cm.dense, norm="log", vmin=4, vmax=22, edgecolors="green", lw=.3)
        ax[i].scatter(D1file.Kpsi[i, 1], D1file.Kphi[i, 1], c=1/D1file.wavenumber, s=6,
                      norm="log", vmin=4, vmax=22, cmap=cm.cm.dense, edgecolors="darkorange", lw=.3)
        ax[i].scatter(D1file.Kpsi[i, 2], D1file.Kphi[i, 2], c=1/D1file.wavenumber, s=6,
                      norm="log", cmap=cm.cm.dense, vmin=4, vmax=22, edgecolors="b", lw=.3)

# Add a colorbar for wavelength
plt.colorbar(ax[-1].collections[0], shrink=.4)

for i in range(3):
    # 1:1 reference line
    ax[i].plot([0, 2], [0, 2], c="k", zorder=0, lw=.5)
    ax[i].set_aspect('equal', adjustable="box")
ax[0].set_xlim(4e-5, 2e0)
ax[0].set_ylim(4e-5, 2e0)

fig.supxlabel(r"K$_\psi$")
fig.supylabel(r"K$_\phi$")
plt.xscale("log")
plt.yscale("log")

plt.legend(fontsize=8, loc="upper right")

plt.savefig("plots/FIGphases.pdf", dpi=500)