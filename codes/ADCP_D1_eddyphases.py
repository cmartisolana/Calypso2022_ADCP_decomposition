# Bin depth criteria of ADCP spectral analysis
# Author: Cristina Martí-Solana
# Description:
#   This script computes mean ADCP spectral densities for different depth layers and eddy phases.
#   For each direction (east/west), it bins spectra by dynamically determined depth indices (e.g., MLD, 200m, 400m)
#   and by eddy phase, then calculates mean spectra and confidence intervals.
#   Results are saved as NetCDF files for reproducibility and publication.

import xarray as xr
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import bootstrap
import statsmodels.stats.api as sms

from functions.open_data import load_section_datasets

import matplotlib.pyplot as plt
import cmocean as cm
plt.rcParams['figure.constrained_layout.use'] = True

DIRECTIONS = ["west", "east"]

mainUCTD = "data/UCTD/"
pathUCTDsections = mainUCTD + "interpolation/"
pathMLD = mainUCTD + "mld/"

def confidence_interval(confidence, data, n, ax):
    """
    Compute confidence intervals for mean spectra.
    """
    mean = np.nanmean(data, axis=ax)
    std = np.nanstd(data, axis=ax)
    if std == 0:
        std = .5 * mean
    lw = mean - confidence * std / np.sqrt(n)
    up = mean + confidence * std / np.sqrt(n)
    return lw, up

for dir in DIRECTIONS:
    # --- Load spectral data for this direction ---
    D0path = f"data/ADCP/D0_files/{dir}.nc"
    D0file = xr.open_dataset(D0path)

    # --- Prepare output file for mean spectra ---
    D1path = f"data/ADCP/D1_files/{dir}.nc"
    D1file = xr.Dataset(
        coords={
            "phase": [0, 1, 2],
            "layer": [0, 1, 2],
            "wavenumber": D0file.wavenumber
        }
    )

    # --- Load section metadata (UCTD and MLD) ---
    pathsUCTD, fileUCTDall = load_section_datasets(path=pathUCTDsections, direction=dir)
    pathsMLD, _ = load_section_datasets(path=pathMLD, direction=dir)

    # --- Calculate mean mixed layer depth and corresponding indices for each transect ---
    MLD = np.empty(len(pathsMLD))
    idxMLD = np.empty(len(pathsMLD))
    for i in range(len(pathsMLD)):
        fileMLD = xr.open_dataset(pathsMLD[i])
        MLD[i] = np.nanmean(fileMLD.depth)
        idxMLD[i] = np.argmin(np.abs(D0file.depth.values - MLD[i]))

    # Indices for 200m and 400m (same for all sections)
    idx200 = np.argmin(np.abs(D0file.depth.values - 200)) * np.ones(len(pathsMLD))
    idx400 = np.argmin(np.abs(D0file.depth.values - 400)) * np.ones(len(pathsMLD))

    # --- Define phase section indices based on direction ---
    if dir == "east":
        idxE = 8
        idxS = 13
    else:
        idxE = 7
        idxS = 12

    # --- Layer and phase indices for binning ---
    # Each entry in layer_idx is an array of indices (per section) for layer boundaries
    layer_idx = [
        np.zeros(len(pathsMLD), dtype=int),    # Surface
        idxMLD.astype(int),                    # MLD
        idx200.astype(int),                    # 200m
        idx400.astype(int)                     # 400m
    ]
    # Phase indices: [start, E, S, end]
    phase_idx = [0, idxE, idxS, len(D0file.section)]

    # --- Compute mean spectra and confidence intervals for each variable ---
    for varname in list(D0file.keys()):
        means_array = np.empty((len(layer_idx) - 1, len(phase_idx) - 1, len(D0file.wavenumber)))
        cilow_array = np.zeros_like(means_array)
        cihigh_array = np.zeros_like(means_array)

        for l in range(len(layer_idx) - 1):  # Loop over layers
            for p in range(len(phase_idx) - 1):  # Loop over phases
                section_means = []
                for s in range(phase_idx[p], phase_idx[p + 1]):  # Loop over sections in this phase
                    start_idx = int(layer_idx[l][s])
                    end_idx = int(layer_idx[l + 1][s])
                    data = D0file[varname][s, start_idx:end_idx, :].values
                    print(f"Processing {varname} for layer {l}, phase {p}, section {s}")
                    # Only compute mean if data is not empty and not all NaN
                    if data.shape[0] > 0 and not np.all(np.isnan(data)):
                        section_mean = np.nanmean(data, axis=0)
                        section_means.append(section_mean)
                # Stack and average across sections, compute confidence intervals
                if section_means and not np.all([np.all(np.isnan(sm)) for sm in section_means]):
                    section_means = np.stack(section_means, axis=0)
                    means_array[l, p, :] = np.nanmean(section_means, axis=0)
                    for k in range(section_means.shape[1]):
                        spatial_data = section_means[:, k]
                        valid_data = spatial_data[~np.isnan(spatial_data)]
                        if valid_data.size < 2:
                            cilow_array[l, p, k] = np.nan
                            cihigh_array[l, p, k] = np.nan
                        else:
                            ci = confidence_interval(0.95, valid_data, len(valid_data), ax=0)
                            cilow_array[l, p, k] = ci[0]
                            cihigh_array[l, p, k] = ci[1]
                else:
                    means_array[l, p, :] = np.nan
                    cilow_array[l, p, :] = np.nan
                    cihigh_array[l, p, :] = np.nan

        # --- Store results in output dataset ---
        D1file[varname] = (["layer", "phase", "wavenumber"], means_array)
        D1file[varname + "_cilow"] = (["layer", "phase", "wavenumber"], cilow_array)
        D1file[varname + "_cihigh"] = (["layer", "phase", "wavenumber"], cihigh_array)

        # Save after each variable for robustness
        D1file.to_netcdf(path=D1path)