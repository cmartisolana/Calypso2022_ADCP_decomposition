# Spectral Helmoltz decomposition for CALYPSO 2022 ADCP data

Spectral decomposition and analysis of ADCP (Acoustic Doppler Current Profiler) data for oceanographic research.

## Overview

This repository contains scripts and data for processing, analyzing, and visualizing ADCP and UCTD (Underway CTD) section data. The workflow includes spectral density computation, Helmholtz decomposition, and the calculation of kinetic and potential energy spectra, with results binned by depth layers and eddy phases.

## Directory Structure

```
codes/
    ADCP_D0_files.py         # Spectral decomposition and Helmholtz analysis
    ADCP_D1_eddyphases.py    # Binning by depth/phase, mean spectra, confidence intervals
    ADCP_D2_FIGphases.py     # Visualization of decomposed spectra
data/
    ADCP/
        ADCP_processed.nc    # Processed ADCP data
        D0_files/            # Spectral results per direction
        D1_files/            # Binned mean spectra per direction
        sections/            # Raw ADCP section files (NetCDF)
    UCTD/
        mld/                 # Mixed Layer Depth data
        sections/            # UCTD section files
plots/                       # Output figures
```

## Main Scripts

- [`codes/ADCP_D0_files.py`](codes/ADCP_D0_files.py): Loads section data, computes spectra, performs Helmholtz decomposition, and saves results to NetCDF.
- [`codes/ADCP_D1_eddyphases.py`](codes/ADCP_D1_eddyphases.py): Bins spectra by dynamically determined depth layers and eddy phases, computes mean spectra and confidence intervals.
- [`codes/ADCP_D2_FIGphases.py`](codes/ADCP_D2_FIGphases.py): Visualizes the relationship between rotational and divergent kinetic energy spectra for different layers and phases.

## Data

- Raw and processed ADCP section data are stored in `data/ADCP/sections/`.
- Processed spectral results are saved in `data/ADCP/D0_files/` and `data/ADCP/D1_files/`.
- UCTD and MLD data are in `data/UCTD/`.

## Requirements

- Python 3.x
- numpy, scipy, xarray, matplotlib, cmocean, scienceplots, statsmodels


## Installation

1. **Install Python dependencies:**
   ```sh
   pip install numpy scipy xarray matplotlib cmocean scienceplots statsmodels
   ```

2. **Install `ADCPySpec` from this repository:**
   Clone or download the `ADCPySpec` module from your repository and ensure it is available in your `PYTHONPATH` or in the `codes/` directory. For example:
   ```sh
   git clone https://github.com/yourusername/ADCPySpec.git
   export PYTHONPATH=$PYTHONPATH:/path/to/ADCPySpec
   ```
   Or copy the `ADCPySpec` directory into your `codes/` folder.

> **Note:** Replace `/path/to/ADCPySpec` with the actual path where you cloned the module.

## Usage

1. **Spectral Decomposition:**  
   Run `ADCP_D0_files.py` to process all ADCP sections and save spectral results.
   ```sh
   python codes/ADCP_D0_files.py
   ```

2. **Binning and Mean Spectra:**  
   Run `ADCP_D1_eddyphases.py` to bin spectra by depth and phase, and compute mean spectra.
   ```sh
   python codes/ADCP_D1_eddyphases.py
   ```

3. **Visualization:**  
   Run `ADCP_D2_FIGphases.py` to generate phase diagrams and other figures.
   ```sh
   python codes/ADCP_D2_FIGphases.py
   ```

## Citation

If you use this code or data in your research, please cite:

Martí-Solana, C. (2025). *Spectral decomposition of ADCP data*. [GitHub Repository](https://github.com/yourusername/ADCP_decomposition)

## License

MIT License

---

For questions or contributions, please open an issue or pull request.
