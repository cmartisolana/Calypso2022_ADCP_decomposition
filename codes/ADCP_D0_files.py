# Spectral decomposition of ADCP data
# Author: Cristina Martí-Solana
# Description:
#   This script processes ADCP and UCTD section data to compute spectral densities,
#   perform Helmholtz decomposition, and calculate potential and kinetic energy spectra.
#   Results are saved as NetCDF files for each direction (east/west).
#   Intended for publication and reproducibility (GitHub-ready).

import numpy as np 
import xarray as xr
import os
import scipy.signal as signal
from scipy.interpolate import interp1d
from scipy.stats import chi2

# ADCPySpec: Custom spectral analysis and Helmholtz decomposition tools
from ADCPySpec.spectrum import SpectrumProcessor
from ADCPySpec.helmholtz import HelmholtzDecomposition

# Internal repository functions
from functions.open_data import load_section_datasets

# Plotting (for diagnostics, not used in output)
import matplotlib.pyplot as plt 
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm
import cmocean as cm
import scienceplots
plt.style.use(['science','no-latex','ieee','high-contrast'])
plt.rcParams['figure.constrained_layout.use'] = True

# --- PARAMETERS ---
DIRECTIONS = ["east","west"]
pathADCP_sections = "data/ADCP/sections"

# Earth's rotation and latitude for Coriolis parameter
omega = 7.2921E-5 # rad/s
lat = 40.76076696859328
f = 2 * omega * np.sin(np.nanmean(lat)*np.pi/180)

# --- MAIN PROCESSING LOOP ---
for dir in DIRECTIONS:
    output_file = f"data/ADCP/D0_files/{dir}.nc"
    paths, all_file = load_section_datasets(direction=dir, path=pathADCP_sections)
    paths = sorted(paths)

    # --- Spectral densities setup ---
    k_grid = np.arange(0, .5, 0.03125)  # Wavenumber grid

    # Preallocate arrays for spectra
    Cuv_grid = np.zeros((len(paths), len(all_file.depth), len(k_grid)))
    Cu_grid = np.zeros_like(Cuv_grid)
    Cv_grid = np.zeros_like(Cuv_grid)

    # --- Compute spectra for each section and depth ---
    for i in range(len(paths)):
        file = xr.open_dataset(paths[i])
        u_track = file.u_track.values
        v_track = file.v_track.values
        distance = file.distance.values

        # Flip distance for east sections to maintain orientation
        if dir == "east":
            distance *= -1
        
        for j in range(len(file.depth)):
            u = u_track[j]
            v = v_track[j]

            # Compute cross-spectra using custom SpectrumProcessor
            spectrum = SpectrumProcessor(distance, u, v, win=None, pad=True)
            k, Cuv, Cu, Cv, _, _, _, _, dof = spectrum.compute_cross_spectrum()

            # Interpolate spectra onto common wavenumber grid
            interCuv = interp1d(k, Cuv, kind='nearest', fill_value="extrapolate")
            Cuv_grid[i, j] = interCuv(k_grid)

            interCu = interp1d(k, Cu, kind='nearest', fill_value="extrapolate")
            Cu_grid[i, j] = interCu(k_grid)

            interCv = interp1d(k, Cv, kind='nearest', fill_value="extrapolate")
            Cv_grid[i, j] = interCv(k_grid)
        
    # --- Mean spectra and confidence intervals ---
    Cuv_mean = np.nanmean(Cuv_grid, axis=0)
    Cu_mean = np.nanmean(Cu_grid, axis=0)
    Cv_mean = np.nanmean(Cv_grid, axis=0)

    def confidence_interval(confidence, data, n, ax):
        """Compute confidence intervals for spectra."""
        mean = np.nanmean(data, axis=ax)
        std = np.nanstd(data, axis=ax)
        lw = mean - confidence * std / np.sqrt(n)
        up = mean + confidence * std / np.sqrt(n)
        return lw, up
    
    Cu_lb, Cu_ub = confidence_interval(.95, Cu_grid, len(k_grid), 0)
    Cv_lb, Cv_ub = confidence_interval(.95, Cv_grid, len(k_grid), 0)

    # --- Helmholtz decomposition ---
    # Preallocate arrays for decomposition results
    Kpsi_grid = np.zeros((len(paths), len(file.depth), len(k_grid)))
    Kphi_grid = np.zeros_like(Cuv_grid)
    Kpsi_iso_grid = np.zeros_like(Cuv_grid)
    Kphi_iso_grid = np.zeros_like(Cuv_grid)

    Cpsi_u = np.zeros_like(Kpsi_grid)
    Cphi_u = np.zeros_like(Kpsi_grid)
    Cpsi_v = np.zeros_like(Kpsi_grid)
    Cphi_v = np.zeros_like(Kpsi_grid)
    Cpsi_uv = np.zeros_like(Kpsi_grid)
    Cphi_uv = np.zeros_like(Kpsi_grid)

    for i in range(len(paths)):
        file = xr.open_dataset(paths[i])
        u_track = file.u_track.values
        v_track = file.v_track.values

        for j in range(len(file.depth)):
            u = u_track[j]
            v = v_track[j]
            theta = np.nanmean(np.arctan(np.abs(v) / np.abs(u)))
            
            decomposition = HelmholtzDecomposition(
                k_grid, Cu_grid[i, j], Cv_grid[i, j], Cuv_grid[i, j], u=u, v=v, theta=theta
            )
            # Isotropic model
            _, _, Kpsi_iso_grid[i, j], Kphi_iso_grid[i, j], _ = decomposition.isotropic_decomposition()
            # Kinetic energy based model for anisotropic flows
            Cpsi_u[i, j], Cphi_u[i, j], Cpsi_v[i, j], Cphi_v[i, j], Cpsi_uv[i, j], Cphi_uv[i, j], Kpsi_grid[i, j], Kphi_grid[i, j], Eu2, Ev2, Euv = decomposition.model3_decomposition()

    # Mean and confidence intervals for Helmholtz components
    Kpsi_iso_mean = np.nanmean(Kpsi_iso_grid, axis=0)
    Kphi_iso_mean = np.nanmean(Kphi_iso_grid, axis=0)
    Kpsi_mean = np.nanmean(Kpsi_grid, axis=0)
    Kphi_mean = np.nanmean(Kphi_grid, axis=0)
    Kpsi_mean[Kpsi_mean == 0] = np.nan
    Kphi_mean[Kphi_mean == 0] = np.nan

    Kpsi_lb, Kpsi_ub = confidence_interval(.95, Kpsi_grid, len(np.where(~np.isnan(Kpsi_grid[:]))), 0)
    Kphi_lb, Kphi_ub = confidence_interval(.95, Kphi_grid, len(np.where(~np.isnan(Kphi_grid[:]))), 0)

    # --- Rossby number as a function of scale ---
    Ro_grid = ((k_grid**3 * .5 * (Cu_grid + Cv_grid) / 1000)**.5) * (1 / f)

    # --- Potential energy spectra from UCTD ---
    # Load UCTD data
    list_paths, fileUCTD_all = load_section_datasets(dir, "data/UCTD/interpolation")
    UCTDdepth = fileUCTD_all.depth.values
    rho0 = np.nanmean(fileUCTD_all.d.values)

    Cb_grid = np.zeros((len(list_paths), len(UCTDdepth), len(k_grid)))
    b_all = []

    for i in range(len(list_paths)):
        fileUCTD_section = xr.open_dataset(sorted(list_paths)[i])
        rho = fileUCTD_section.d.values
        g = 9.81  # Gravitational acceleration (m/s²)

        # Buoyancy
        b_track = -g * (rho - rho0) / rho0
        b_track[np.isnan(rho)] = 0.
        b_all.extend(b_track.flatten())

        # Brunt-Väisäla frequency (N^2)
        depthgradient = np.gradient(rho, -4, axis=0)
        N_track = ((-g / rho0) * depthgradient)**(1/2)

        for j in range(len(UCTDdepth)):
            b = b_track[j]
            N = np.nanmean(N_track[j])

            spectrum = SpectrumProcessor(distance, b, win=None, pad=True)
            k, Cb = spectrum.compute_one_signal_spectrum()
            Cb = (Cb / N)
            interCb = interp1d(k, Cb, kind='nearest', fill_value="extrapolate")
            Cb_grid[i, j] = interCb(k_grid)
    
    # Interpolate Cb_grid to match the depth of the ADCP data
    Cb_grid_interp = np.zeros((len(paths), len(file.depth), len(k_grid)))
    for i in range(len(paths)):
        for k in range(len(k_grid)):
            interp_func = interp1d(UCTDdepth, Cb_grid[i, :, k], kind='linear', fill_value="extrapolate")
            Cb_grid_interp[i, :, k] = interp_func(file.depth)

    # --- Wave-vortex decomposition ---
    E = .5 * (Cu_grid + Cv_grid + Cb_grid_interp)
    Ew = 2 * Kphi_grid
    Ev = E - Ew

    # --- Save all results to NetCDF ---
    file_spectrum = xr.Dataset(
        data_vars={
            "Cpsi_u": (["section", "depth", "wavenumber"], Cpsi_u),
            "Cphi_v": (["section", "depth", "wavenumber"], Cphi_v),
            "Kphi_iso": (["section", "depth", "wavenumber"], Kphi_iso_grid),
            "Kphi": (["section", "depth", "wavenumber"], Kphi_grid),
            "Kpsi_iso": (["section", "depth", "wavenumber"], Kpsi_iso_grid),
            "Kpsi": (["section", "depth", "wavenumber"], Kpsi_grid),
            "KE": (["section", "depth", "wavenumber"], (Cu_grid + Cv_grid) / 2),
            "PE": (["section", "depthUCTD", "wavenumber"], Cb_grid / 2),
            "E": (["section", "depth", "wavenumber"], E),
            "Ew": (["section", "depth", "wavenumber"], Ew),
            "Ev": (["section", "depth", "wavenumber"], Ev),
            "Ro": (["section", "depth", "wavenumber"], Ro_grid)
        },
        coords={
            "section": range(len(paths)),
            "depth": file.depth,
            "depthUCTD": UCTDdepth,
            "wavenumber": k_grid
        }
    )

    file_spectrum.to_netcdf(path=output_file)