# CP-Analyzer: Constant-Potential Free Energy Analysis for Electrocatalysis

**CP-Analyzer** is an open-source Python package for post-processing VASP calculations to perform constant-potential electrocatalytic free energy analysis.It enables systematic evaluation of adsorption energies, reaction pathways, and catalytic activity across a wide range of electrode potentials and pH conditions, based directly on VASP outputs.Unlike traditional constant-charge approaches, CP-Analyzer implements a quadratic fitting strategy derived from constant-potential DFT data, allowing reliable free-energy prediction under realistic electrochemical conditions. The package provides a complete workflow from raw VASP data extraction to microkinetic analysis and figure generation.

## Key Features

- Automated extraction of data from VASP outputs
- Constant-potential free energy evaluation via quadratic fitting
- Analysis across the full pH–potential space
- Automatic conversion between SHE and RHE scales
- Microkinetic analysis tools
- Publication-quality plotting utilities
- Modular task-based architecture for flexible execution

## Contents

- Citation
- Installation
- Usage
- Core Code Structure
- Examples
- Authors
- License

## How to cite

If CP-Analyzer contributes to your research, please cite:

CP-Analyzer: A Comprehensive Platform for Constant-Potential Free Energy Calculations and Analysis in Electrocatalysis

## Prerequisites

### Operating System

- Linux / macOS / Windows (64-bit)

### Python

- Python ≥ 3.8

### Dependencies (tested versions)

CP-Analyzer has been tested with the following package versions:

- matplotlib == 3.10.9
- numpy == 2.4.4
- pandas == 3.0.2
- scipy == 1.17.1
- setuptools == 82.0.1
- Unidecode == 1.4.0

> Other nearby versions may also work, but the above configuration is verified to be fully compatible.

## Installation (Recommended via Conda)

The easiest way to install CP-Analyzer is by using Conda to create a clean Python environment.

### 1. Create and activate a new environment

```bash
conda create -n cp-analyzer python=3.10 -c conda-forge
conda activate cp-analyzer
```

### 2. Clone the repository

```bash
git clone https://github.com/lxf-gzu/CP-Analyzer.git
cd CP-Analyzer
```

### 3. Install CP-Analyzer and dependencies

```python
pip install .
```

> This command automatically installs all required dependencies and registers the cp-analyzer command.

#### Test Installation

```python
cp-analyzer -h
```

> If the help message is displayed, the installation is successful.

#### Deactivate Environment

```bash
conda deactivate
```

## Core Package Structure

The cp_analyzer package contains the core modules of the platform:

- **utils.py** — Utility functions
- **data_processing.py** — VASP data extraction and preprocessing
- **calculations.py** — Constant-potential free energy calculations
- **plotting.py** — Figure generation
- **tasks.py** — Analysis task implementations
- **main.py** — Command-line entry point
- **config_default.py** — Default configuration template (automatically copied to working directory as config.py when first executed)

## Tips for Configuration (config.py)
Some example datasets are provided in the CP-Analyzer/examples directory, including slab.xlsx, ooh.xlsx, o.xlsx, oh.xlsx, and CO.xlsx, together with representative VASP input and output files (e.g., INCAR, POSCAR, KPOINTS, OUTCAR, and LOCPOT) for demonstration and testing purposes. These examples cover typical adsorption and reaction intermediate systems and can be directly used to validate the workflow.

In addition, CP-Analyzer includes two built-in datasets, literature and literature_params, which can be used for quick testing without requiring any external files.

### Data Source Settings

```python
# General output settings
DEFAULT_OUTPUT_DIR = "FreeEnergy_Output_file"  # Default directory for all output files
# Input file sources
DEFAULT_DATA_SOURCE = "file" # Source: "file", "literature", or "literature_params"
DEFAULT_FILE_PATHS = ['slab.xlsx', 'ooh.xlsx', 'o.xlsx','oh.xlsx']
# Input Excel files
```

### Zero-Point Energy (ZPE) Corrections

```python
# Zero-point energy corrections for literature and file-based data
DEFAULT_E_ZPE_ST_LITERATURE = [0, 0.3257, 0.015107, 0.291982] #(Fig.3 c)https://doi.org/10.1021/jacs.2c08743 with zero point energy correction
#DEFAULT_E_ZPE_ST_LITERATURE = [0, 0, 0, 0] # (Fig.3b d)https://doi.org/10.1021/jacs.2c08743 has no zero point energy correction
#DEFAULT_E_ZPE_ST_LITERATURE = [0, 0, 0, 0]  # For literature data,slab->ooh->o->oh
DEFAULT_E_ZPE_ST_FILE = [0, 0.334185, -0.00551, 0.308282]  # For file-based data,slab->ooh->o->oh
```

The zero-point energy (ZPE) corrections defined in DEFAULT_E_ZPE_ST_LITERATURE and DEFAULT_E_ZPE_ST_FILE are derived from our dataset. Some examples use experimental data from Hu et al.1, however, since that work does not provide corresponding free energy corrections, we apply the ZPE corrections from our own dataset. As a result, certain fitted results may show slight deviations from those reported in that work, but remain within a reasonable range. Users can modify these values as needed.

The electrode potential scale can be adjusted, the onset potential at various pH values can be computed, and calculations involving O₂ intermediates can be performed in this section:

### Onset Potential Settings

```python
DEFAULT_RUN_ONSET_ANALYSIS = True     # Enable/disable Task 4 execution
DEFAULT_VOLTAGE_SCALE = "RHE"         # Voltage scale: "RHE" or "SHE" (must be RHE for onset calc)
DEFAULT_ONSET_PH_VALUES = [1, 13]     # pH values for onset potential calculations
DEFAULT_INCLUDE_O2 = False            # Include O2 in reaction intermediates
```

### Example:

```python
DEFAULT_RUN_ONSET_ANALYSIS = True 
DEFAULT_VOLTAGE_SCALE = "RHE" 
DEFAULT_ONSET_PH_VALUES = [1, 7，13] 
DEFAULT_INCLUDE_O2 = True 
```

Users can make modifications here if they need to use the free energy profiles at different pH values.

### Single pH Free Energy Analysis

```python
DEFAULT_RUN_SINGLE_PH = True          # Enable/disable Task 5 execution
DEFAULT_ANALYSIS_PH_VALUES = [1]      # pH values for single point analysis
# Experimental RHE voltage ranges for OER/ORR
EXPERIMENTAL_RHE_RANGES = {
    "ORR": (-3, 2),     # Voltage range for Oxygen Reduction Reaction
    "OER": (1.2, 2.6)   # Voltage range for Oxygen Evolution Reaction
}
```

> **Note**: DEFAULT_ANALYSIS_PH_VALUES = [7] indicates analysis at neutral pH, which can be customized by users. EXPERIMENTAL_RHE_RANGES represents the experimental potential ranges for ORR and OER.
Main Program (main.py)

### Task Dependencies

```python
def main():
    TASK_DEPENDENCIES = {
        2: [1],
        3: [1],
        4: [1],
        5: [1],
        6: [1],
        7: [1],
        8: [1],
        9: [1],
        10: [1],
        13: [1],
        12: [11],
    }
```

> **Note**:Tasks 2–10 and 13 depend on Task 1 (fitting parameters) ,and Task 12 depends on Task 11 (VASP directory setup).

## Usage

### Test with Built-in Data:

```python
cp-analyzer  --data-source literature --task 1 2 4
```

### Use Custom VASP Data:

```python
cp-analyzer  --data-source file \ 
--file-paths slab.xlsx ooh.xlsx o.xlsx oh.xlsx \ 
--task 1 2 4 5 6 7 \ 
--voltage-scale RHE \ 
--output-dir my_analysis
```

### Run Full Workflow:

```python
cp-analyzer  --data-source file \ 
--file-paths slab.xlsx ooh.xlsx o.xlsx oh.xlsx \ 
--task 1-15 \ 
--voltage-scale RHE \ 
--onset-ph-values 0 7 14 \ 
--analysis-ph-values 1 13 \ 
--output-dir full_analysis
```

## Detailed Features

To fully use this workflow, please read this section carefully. Users should first run Task 11 to create VASP calculation directories with different charge states, after VASP calculations finished, run Task 12 to extract VASP data. The workflow is as follows:

### Task 11: VASP Directory Creation

**Function**: Batch creation of VASP calculation directories with different charge states. To use this script, place all required VASP input files in the current directory. In the INCAR file, the NELECT tag must be commented out and written in the form # NELECT = xxx, where xxx is the total number of electrons of the neutral system obtained from the OUTCAR (you can use grep “NELECT” OUTCAR). This commented line is required because the script reads it to determine the neutral reference electron number and automatically adjusts NELECT for different charge states when generating calculation directories, as shown below:

```bash
dir/ 
└── POSCAR/ 
├── INCAR 
├── POTCAR/ 
├── KPOINTS/ 
├── jobs.vasp/ #提交任务脚本
```

**Commands**:

```python
# create directories only
cp-analyzer  --task 11
# create directories and submit jobs automatically
cp-analyzer  --task 11 --submit-vasp-jobs
# custom input files
cp-analyzer  --task 11 \
--vasp-input-files INCAR KPOINTS POSCAR POTCAR submit.sh \
--vasp-submit-script submit.sh
```

**Output**:

```bash
FreeEnergy_Output_file/
└── TASK_11_vasp_directories/
├── e0/ # neutral system
├── e-0.5/ # -0.5e charge
├── e-1/
├── e-1.5/
├── e-2/
├── e+0.5/ # +0.5e charge
├── e+1/
├── e+1.5/
└── e+2/
```

### Task 12: VASP Data Extraction

**Function**: Extract electrochemical data from VASP output files (LOCPOT, OUTCAR). This task depends on the LOCPOT and OUTCAR files generated in Task 11. Run the following command:

```python
cp-analyzer  --task 12
```

**Output**:

```bash
FreeEnergy_Output_file/
└── data_for_different_charges.xlsx
```

The file data_for_different_charges.xlsx contains the following five types of data:

- **E_fermi**: Fermi level (from OUTCAR)
- **top_E_Vacums**: vacuum level at the top (from planar-averaged LOCPOT)
- **E_potentials**: average electrostatic potential (from LOCPOT)
- **E_DFTs**: total DFT energy TOTEN (from OUTCAR)
- **Charges**: system charge (corresponding to directory charge increment)

### Task 1: Parameter Fitting

**Function**: Fit quadratic parameters (C, U_PZC, E0, R²) from constant-potential DFT data Multiple data sources are supported. The first is the dataset from Hu et al. [1], including "literature" and "literature_params", which are built into config.py. The second is the dataset generated from VASP calculations, "file" (slab.xlsx, ooh.xlsx, o.xlsx, oh.xlsx). Users should place these files in the root directory and set the data source in config.py:

**Set Data Source**

```python
# Set Data Source
DEFAULT_DATA_SOURCE = "literature "# Source: "file", "literature", or "literature_params"
```

**Methodology**:The code reads electrochemical data from Excel files (or built-in literature data), calculates free energies versus SHE potential using E_free_vs_SHE, and performs bounded nonlinear least-squares fitting to the quadratic function:

$$ E_i(U_{\rm SHE}) = -\frac{1}{2} C_i (U_{\rm SHE} - U_{\rm PZC,i})^2 + E_{0,i} $$

where idenotes slab, *OOH, *O, *OH, etc. The fitting uses robust initial guesses and bounds to ensure physical results (C > 0).

**Supported Data Sources**:
- "file" (recommended for user VASP data)
- "literature" / "literature_params" (for quick testing)

**Command**:

```python
cp-analyzer  --task 1
```

**Main Outputs**:
- Summary table of fitting parameters (C, U_PZC, E0, R²)
- Individual fitting information for each species

> **Note**: All subsequent tasks (except Task 11, 12, 14) depend on the successful completion of Task 1.

### Task 2: Free Energy–Potential Curve

**Function**: Visualizes the fitted quadratic free energy curves as a function of electrode potential versus SHE for all surface states, allowing intuitive understanding of potential-dependent thermodynamics.

**Methodology**:Using the parameters obtained in Task 1, the code generates continuous energy curves over a wide potential range and overlays original data points.

**Command**:

```python
cp-analyzer  --task 1 2
```

**Main Outputs**:
- free_energy_vs_she.png (publication-quality plot)
- Raw and fitted data for each species (stored in free_energy_fits/ folder)

### Task 3: Q–pH Analysis

**Function**: Task 3 analyzes the relationship between the charge transfer (denoted as dQ) of the substrate and adsorbed intermediates and the solution pH. The main goal is to reveal the coupling between proton transfer (pH) and electron transfer at the electrocatalytic interface under constant-potential conditions. By constructing continuous dQ–pH curves, this module helps understand how the charge state of different reaction intermediates evolves with pH, providing mechanistic insight into proton-coupled electron transfer (PCET) processes.

**Methodology**: For each input Excel file corresponding to a specific surface state (slab, *OOH, *O, *OH, etc.), the code performs the following steps:

1. Reads the electrochemical data from the first five columns and computes the electrode potential versus the Standard Hydrogen Electrode ($E_U$ vs SHE ) using the function E_free_vs_SHE.

2. Calculates the corresponding pH value for each data point using the Nernst relation (with $U_{RHE}=0$ V as reference): $$\mathrm{pH} = \frac{-E_{U}\ \rm vs\ SHE}{0.0592}$$ where 0.0592 V is the Nernst factor (2.303 RT/F) at 298 K.

3. Uses the user-provided dQ value (sixth column) for that data point.

4. Collects pairs of (𝑝𝐻, 𝑑𝑄) for each surface state and performs a third-order polynomial fitting: $$\mathrm{d}Q(\mathrm{pH}) = a \cdot \mathrm{pH}^3 + b \cdot \mathrm{pH}^2 + c \cdot \mathrm{pH} + d$$

5. Generates fitting curves, predicted dQ values over pH = 0–14, and saves the corresponding figures and data tables.

**Input File Preparation** (Critical)

To execute Task 3, each Excel file (slab.xlsx, ooh.xlsx, o.xlsx, oh.xlsx, etc.) must contain six columns with the following structure:

- Columns 1–5: Standard electrochemical data required by the platform (E_fermi, top_E_Vacuum, E_potential, E_DFT, Charge).
- Column 6 (last column): dQ — the charge transfer amount.

**Definition of dQ**:
- dQ represents the net charge difference associated with the adsorption of the intermediate relative to the bare slab in the same charge state.
- For the bare slab.xlsx: dQ = 0 (no adsorbate, no charge transfer).
- For adsorbate systems (*OOH, *O, *OH, etc.): dQ should be the charge variation caused by the presence of the adsorbed intermediate (commonly obtained from Bader charge analysis or charge density integration)

Each row in the Excel file corresponds to one DFT calculation with a specific charge state (different NELECT). The code automatically pairs the computed pH (from the potential) with the user-supplied dQ to construct the pH-dependent charge response.

**Command**:

```python
cp-analyzer  --task 1 3
```

**Output**
- Q_vs_pH_fit.png: Fitted dQ–pH curves for all species.
- Q_vs_pH_range.png: Predicted dQ values across pH 0–14.
- Text files containing fitting parameters and raw/predicted data.

> **Note**: If the sixth column is missing or contains no valid numeric values, Task 3 will be automatically skipped with a warning message.

### Task 4: Onset Potential Analysis

**Function**: Calculates the onset potentials and identifies the rate-determining steps (RDS) for ORR and OER at different pH values. It accounts for the fact that the RDS may shift with applied potential.

**Methodology**: The code first determines the initial RDS at U = 0 VRHE. It then uses a self-consistent multi-stage algorithm (coarse grid search + bisection method + local optimization) to find the potential Uwhere the maximum elementary step free energy satisfies max⁡(ΔGi(U))≈0. Full step energies are recomputed at the solved onset potential to confirm the actual RDS. The details see the CP-Analyzer related paper.

**Command**:

```python
cp-analyzer  --task 1 4
```

**Main Outputs**:

- onset_potentials.csv and onset_potentials.json
- Detailed summary table with pH, RDS, and onset potential (vs RHE)

> **Note**: Requires DEFAULT_VOLTAGE_SCALE = "RHE".

### Task 5: Single pH Analysis

**Function**: Computes detailed potential-dependent free energy profiles (intermediates and four elementary steps) for ORR and OER at user-specified pH values.

**Methodology**: For each pH, the code evaluates adsorption free energies using the quadratic model across an experimental potential range, then calculates step free energies for both ORR and OER.

**Command**:

```python
cp-analyzer  --task 1 5
```

**Main Outputs**:
- CSV files containing intermediates and step energies vs potential
- Plots of intermediate energies and step energies at each pH

### Task 6: Heatmap Analysis

**Function**: Creates two-dimensional heatmaps showing how free energies of intermediates and elementary steps vary across the full pH–potential space.

**Methodology**: Calculates properties on a dense grid (default 150×150) of pH and potential values using the quadratic model, then generates heatmaps for visualization.

**Command**:

```python
cp-analyzer  --task 1 6
```

**Main Outputs**:

- Heatmap figures for ORR and OER (intermediates + steps) in the heatmaps/ folder

### Task 7: Reaction Pathway Diagram

**Function**: Plots free energy diagrams of the ORR and OER reaction pathways under multiple pH and electrode potential conditions, with automatic identification of the rate-determining step (RDS) for each condition.

**Methodology**: The code uses the quadratic model parameters from Task 1 to compute the free energies of intermediates and elementary steps, then constructs cumulative free energy profiles along the reaction pathway.

**Important Note (Recommended Usage)**: It is highly recommended to first run Task 4 (Onset Potential Analysis) to obtain accurate onset potentials, and then update the conditions in config.py to include these onset potentials for more meaningful pathway diagrams.

```python
# Conditions for pathway diagrams in config.py
DEFAULT_ORR_STEP_CONDITIONS = [
    {'label': 'pH=1, U=0 RHE (ORR)', 'pH': 1, 'U': 0},
    {'label': 'pH=13, U=0 RHE (ORR)', 'pH': 13, 'U': 0},
    
    # Use onset potentials obtained from Task 4
    {'label': 'pH=1, U=0.5077 RHE (ORR Onset)', 'pH': 1, 'U': 0.5077},
    {'label': 'pH=13, U=0.3344 RHE (ORR Onset)', 'pH': 13, 'U': 0.3344}
]

DEFAULT_OER_STEP_CONDITIONS = [
    {'label': 'pH=1, U=0 RHE (OER)', 'pH': 1, 'U': 0},
    {'label': 'pH=13, U=0 RHE (OER)', 'pH': 13, 'U': 0},
    
    # Use onset potentials obtained from Task 4
    {'label': 'pH=1, U=1.7330 RHE (OER Onset)', 'pH': 1, 'U': 1.7330},
    {'label': 'pH=13, U=1.8771 RHE (OER Onset)', 'pH': 13, 'U': 1.8771}
] 
```

**Command**:

```python
cp-analyzer  --task 1 7
```

**Main Outputs**:

- Reaction pathway diagrams (PDF/PNG) in the reaction_pathway_ORR/ and reaction_pathway_OER/ folders
- Cumulative free energy profiles with highlighted RDS for each condition

> **Tip**: Running Task 4 first and then updating the U values in the conditions above will generate diagrams at both equilibrium (U=0) and onset potentials, providing deeper mechanistic insight.

### Task 8: Microkinetic Parameter Export

**Function**: Exports the free energy changes (ΔG) of the elementary steps for ORR and OER at specified pH and electrode potential conditions. These parameters are primarily used as input for Task 13 (Microkinetic Simulation) or other external microkinetic modeling software.

**Methodology**: Using the quadratic model parameters from Task 1, the code calculates the reaction free energies (ΔG1 to ΔG4) for ORR and OER at user-defined pH values and potentials (default: pH = 1 and 13, U = 0 VRHE). The results are saved in a simple text format.

```python
# MICROKINETICS Parameters for ORR (RHE voltage)
# Format: pH, Potential vs RHE (V), dG1, dG2, dG3, dG4
1.0000, 0.0000, -0.618389, -1.552131, -1.701575, -1.047905
13.0000, 0.0000, -0.989798, -1.930957, -1.295692, -0.703553
```

**Key Output Example** (microkinetics_ORR.txt):

**Important Connection with Task 13**: The four ΔGvalues (dG1 to dG4) for ORR generated by Task 8 are directly used to update the "delta_Gi" parameter in MICROKINETICS_CONFIG of config.py (the last four values, excluding the first two zeros).

**Command**:

```python
cp-analyzer  --task 1 8
```

**Main Outputs**:

- microkinetics_ORR.txt
- microkinetics_OER.txt

> **Note**: It is recommended to run Task 8 after Task 1, then copy the generated ΔG values into MICROKINETICS_CONFIG["delta_Gi"] before running Task 13 for consistent microkinetic simulations.

### Task 9: Charge–Potential Analysis

**Function**: Analyzes the linear relationship between system charge (Q) and electrode potential (U), extracts interfacial capacitance, and provides practical recommendations for setting NELECT in VASP calculations.

**Methodology**: Performs linear regression $Q=C⋅U+Q_0$ on the charge-potential data and evaluates charge at key potentials.

**Command**:

```python
cp-analyzer  --task 1 9
```

**Main Outputs**:

- Linear fitting results and capacitance values
- Text files with NELECT adjustment suggestions at typical potentials (-1.23 V, 0 V, 1.23 V)

### Task 10: General Adsorption Systems

**Function**: Extends the quadratic fitting framework to arbitrary adsorbates (e.g., *CO, *H, etc.), allowing flexible analysis of potential-dependent adsorption energies.

**Input Preparation**: Place corresponding Excel files (e.g., CO.xlsx) in the root directory and update ADS_FILE_PATHS and ADS_SYSTEMS in config.py.

**Command**:

```python
cp-analyzer  --task 1 10
```

**Main Outputs**:

- Fitting parameters and adsorption energy vs potential plots

### Task 13: Microkinetic Simulation

**Function**: Performs full microkinetic modeling for ORR and OER based on the free energies obtained from the quadratic model. It solves the coupled differential equations under the steady-state approximation to simulate surface coverages of intermediates and macroscopic current density as a function of electrode potential, providing a bridge from atomic-scale energetics to observable catalytic performance.

**Methodology**: The module uses mean-field microkinetic modeling. It constructs a reaction network consisting of 6 elementary steps (including adsorption, desorption, and electrochemical proton-coupled electron transfer steps). The rate constants are calculated using transition state theory and the Butler–Volmer formalism. The system of ordinary differential equations for intermediate coverages is solved numerically at steady state over a range of potentials. The details please see the CP-Analyzer related paper.

**Key Parameters** (configured in config.py under MICROKINETICS_CONFIG):

```python
MICROKINETICS_CONFIG = {
    "acc": 50,                          # Numerical precision for symbolic solver
    "x_H2O": 1.0,                       # Mole fraction of water
    "x_O2": 2.34e-5,                    # Mole fraction of O2
    "T": 298.0,                         # Temperature (K)
    "beta_i": [0, 0, 0.5, 0.5, 0.5, 0.5],   # Electron transfer coefficients for each step
    "delta_Ga0i": [0.0, 0.0, 0.260, 0.260, 0.260, 0.260],  # Intrinsic activation barriers (eV)
    "delta_Gi": [0.0, 0.0, -0.618389, -1.552131, -1.701575, -1.047905],  # Reaction free energies (eV) — usually updated from Task 8
    "Ai": [8.0e5, 1.0e8, 1.0e9, 1.0e9, 1.0e9, 1.0e9],   # Pre-exponential factors (s⁻¹)
    "current_max": 15,                  # Maximum current density for plot truncation (mA/cm²)
    "U_range_start": 0.0001,            # Start potential (V vs RHE)
    "U_range_stop": 2.0,                # End potential (V vs RHE)
    "U_range_step": 0.01                # Potential step size (V)
}
```

**Important Workflow**:

1. Run Task 1 to obtain quadratic fitting parameters.
2. Run Task 8 to generate ΔGvalues at desired pH.
3. Copy the four ΔG values from Task 8 output into MICROKINETICS_CONFIG["delta_Gi"] (the last four elements).
4. Run Task 13.

**Command**:

```python
cp-analyzer  --task 1 13
```

**Main Outputs** (saved in the root output directory or Task13 folder):

- Polarization curves (ORR_OER_Polarization_Curve.png)
- Surface coverage plots (Intermediate_Coverages.png, O2_dl_Coverage.png)
- TOF (Turnover Frequency) plots
- Text files: U_coverage_current_data.txt, Key_potentials.txt

> **Note**: Task 13 is computationally more intensive than other tasks. The accuracy of the simulation strongly depends on the quality of the free energies (delta_Gi) provided from Task 8 and the chosen activation barriers (delta_Ga0i).

### Task 14: Constant-Charge Reaction Diagram

**Function**: Generates reaction free energy diagrams (stepwise energy profiles) for ORR and OER directly from constant-charge VASP calculation results, without applying quadratic fitting. This task allows users to quickly visualize and compare free energy landscapes at specific pH and electrode potential (U) values under the traditional constant-charge approximation.

**Methodology**: The module uses the raw DFT energies defined in OER_ORR_CONFIG (in config.py), applies thermodynamic corrections (ZPE, entropy, pH, and potential), and constructs the four-step free energy diagram for the chosen reaction (ORR or OER). It supports flexible specification of pH and applied potential U.

**Configuration** (in config.py):

```python
DEFAULT_RUN_OER_ORR = True      # Enable/disable Task 14 execution

# Default configurations for OER/ORR pathway analysis
OER_ORR_CONFIG = {
    "pH": 0.0,                  # Default pH value
    "U": 0.0,                   # Default potential (V vs. SHE)
    "reaction_type": "ORR",     # Default reaction: "ORR" or "OER"
    # VASP energies (example values - should be replaced with actual calculations)
    "vasp_energies": {
        "E_slab": -100.0,       # Total energy of substrate (eV)
        "E_OOH": -120.0,        # Total energy of *OOH (eV)
        "E_O": -110.0,          # Total energy of *O (eV)
        "E_OH": -105.0,         # Total energy of *OH (eV)
        "E_H2O": -14.22,        # Total energy of H2O (eV)
        "E_O2": -9.86,          # Total energy of O2 (eV)
        "E_H2": -6.77           # Total energy of H2 (eV)
    },
    # Energy corrections
    "corrections": {
        "delta_G_water": 4.92,  # H2O formation free energy (eV)
        "ZPE_TDS": {            # Zero-point energy and thermal corrections (eV)
            "*OOH": 0.2,
            "*O": 0.1,
            "*OH": 0.15
        }
    },
    "U_range": np.linspace(0.0, 1.5, 100),  # Voltage range (for reference)
    "equilibrium_potential": 1.23           # Standard potential for ORR/OER
}
```

**Command**:

```python
# Basic usage
cp-analyzer  --task 14

# Example: ORR at pH=1, U=0.8 V vs SHE
cp-analyzer  --task 14 --pH 1 --U 0.8 --reaction-type ORR
```

**Main Outputs**:

- Free energy step diagrams (PDF/PNG)
- CSV files containing intermediate energies and elementary step free energies (ΔG1, ΔG2, ΔG3, ΔG4)
- RDS (rate-determining step) information for the given condition

> **Note**:
- Unlike most other tasks, Task 14 does not require quadratic fitting (Task 1). It directly uses the energies provided in OER_ORR_CONFIG.
- This task is particularly useful for rapid comparison between constant-charge and constant-potential (quadratic model) approaches.

### Task 15: Activity Volcano Plot

**Function**: Constructs two-dimensional activity volcano plots for ORR and OER based on the Sabatier principle. It uses adsorption free energies of key intermediates (*OH, *O, *OOH) as descriptors to predict catalytic activity (overpotential) across a wide range of materials, helping users identify optimal catalysts and understand scaling relations.

**Methodology**: The module performs the following steps:
1. Uses a built-in database of adsorption free energies for various metals (METAL_ADSORPTION_ENERGIES in config.py).
2. Applies linear scaling relations (especially between *OOH and *OH) derived from the data.
3. Calculates the free energy changes of all four elementary steps for ORR or OER on a dense 2D grid of descriptor values ($ΔG_{OH^∗}$ and $ΔG_{O^∗}−ΔG_{OH^∗}$).
4. Determines the theoretical overpotential ($η$) as the minimum potential required to make all steps downhill (i.e., max⁡($ΔG_i$)−1.23 V).
5. Generates contour plots with metal data points, boundary lines separating different rate-determining steps, and exports all underlying data.

```python
# Metal adsorption energies: (ΔG_OH, ΔG_O, ΔG_OOH)
METAL_ADSORPTION_ENERGIES = {
    "Ti": (-0.93, -0.74, 2.5),
    "Pt": (0.19, 2.33, 3.56),
    # ... more metals
}

VOLCANO_PLOT_RANGE = (-1.5, 3.0, -0.5, 3.5)   # (x_min, x_max, y_min, y_max)
OVERPTENTIAL_COLOR_RANGE = (0.32, 3.0)  # Color scale for overpotential
```

**Configuration** (in config.py):

**Command**:

```python
cp-analyzer  --task 15
```

**Main Outputs** (saved in TASK_15_oer_orr_volcano/2D_Volcano/):
- ORR_activity_volcano.png and OER_activity_volcano.png
- Grid data (*_volcano_grid_data.csv)
- Metal points data (*_volcano_metal_points.csv)
- Boundary line data separating different RDS regions
- Linear scaling fit parameters

## Data

All example datasets in this package are available in the root directory.

## Authors

This software package is mainly developed and validated by Prof. Xuefei Liu.

## License

CP-Analyzer is released under the MIT License.

## Reference

(1) Hu, X.; Chen, S.; Chen, L.; Tian, Y.; Yao, S.; Lu, Z.; Zhang, X.; Zhou, Z. What is the real origin of the activity of Fe–N–C electrocatalysts in the O2 reduction reaction? Critical roles of coordinating pyrrolic N and axially adsorbing species. Journal of the American Chemical Society 2022, 144 (39), 18144-18152.
