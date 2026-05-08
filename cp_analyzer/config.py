#########################################################
"""
Electrocatalysis Simulation Toolkit
Centralized configuration for ORR/OER simulations
Core Functionalities:
1. Constant Potential Calculations
2. Microkinetic Modeling
3. Routine Catalysis Data Processing
"""
#########################################################

import numpy as np

############################
# TASK-INDEPENDENT SETTINGS
############################

# General output settings
DEFAULT_OUTPUT_DIR = "FreeEnergy_Output_file"  # Default directory for all output files
# Input file sources
DEFAULT_DATA_SOURCE = "file"      # Source: "file", "literature", or "literature_params"
# literature doi.org/10.1021/jacs.2c08743 Figure 5, literature_params for figure 3
DEFAULT_FILE_PATHS = ['slab.xlsx', 'ooh.xlsx', 'o.xlsx', 'oh.xlsx']  # Input Excel files

###########################
# ZPE CORRECTIONS
# Used by multiple tasks (4-8, 14)
###########################

# Zero-point energy corrections for literature and file-based data
DEFAULT_E_ZPE_ST_LITERATURE = [0, 0.3257, 0.015107, 0.291982] #(Fig.3 c)https://doi.org/10.1021/jacs.2c08743 with zero point energy correction
#DEFAULT_E_ZPE_ST_LITERATURE = [0, 0, 0, 0] # (Fig.3 b d)https://doi.org/10.1021/jacs.2c08743 has no zero point energy correction
#DEFAULT_E_ZPE_ST_LITERATURE = [0, 0, 0, 0]  # For literature data,slab->ooh->o->oh
DEFAULT_E_ZPE_ST_FILE = [0, 0.334185, -0.00551, 0.308282]  # For file-based data,slab->ooh->o->oh

# Thermodynamic constants for free energy calculations
G_O2 = -9.903      # Free energy of O2 gas (eV) at standard conditions
G_H2 = -6.811      # Free energy of H2 gas (eV) at standard conditions
G_H2O = -14.223    # Free energy of H2O liquid (eV) at standard conditions
T = 298.15         # Temperature in Kelvin
k = 8.617385e-5    # Boltzmann constant in eV/K
delta_G_water = 4.92  # Free energy of water formation (eV)


########################
# TASK 1-2 CONFIGURATION
# Used for parameter fitting and free energy plots
########################

# Literature data for electrochemical systems
# Dict keys: system names, values: dict with 'U' (potential vs SHE) and 'E' (energy in eV) arrays
# doi.org/10.1021/jacs.2c08743 Figure 5
LITERATURE_DATA = {
    "slab": {
        "U": np.array([-3.01, -2.12, -1.42, -0.46, 0.05, 0.55, 1.01, 1.45, 1.86]),
        "E": np.array([-227.93, -226.38, -225.56, -224.86, -224.75, -224.86, -225.23, -225.77, -226.54])
    },
    "*OOH": {
        "U": np.array([-2.01, -1.27, -1.05, -0.35, 0.14, 0.63, 1.10, 1.42, 1.88]),
        "E": np.array([-241.94, -241.15, -240.36, -239.77, -239.69, -239.81, -240.15, -240.46, -241.40])
    },
    "*O": {
        "U": np.array([-2.09, -1.47, -0.54, -0.09, 0.40, 0.86, 1.24, 1.72, 2.10]),
        "E": np.array([-233.22, -232.12, -231.03, -230.69, -230.61, -230.69, -230.98, -231.60, -232.31])
    },
    "*OH": {
        "U": np.array([-2.21, -1.78, -1.10, -0.40, 0.18, 0.73, 1.10, 1.49, 1.94]),
        "E": np.array([-237.84, -236.68, -235.83, -235.27, -235.16, -235.11, -235.65, -236.13, -236.93])
    }
}

# Pre-fitted parameters for literature systems
# Used when data_source = "literature_params"
# doi.org/10.1021/jacs.2c08743 Figure 3
LITERATURE_PARAMS = {
    "slab": {"C": 1.074, "U_pzc": -0.773, "E0": -281.089},
    "*OOH": {"C": 1.031, "U_pzc": -0.587, "E0": -296.503},
    "*O": {"C": 0.959, "U_pzc": -0.470, "E0": -287.539},
    "*OH": {"C": 0.978, "U_pzc": -0.649, "E0": -291.903},
}

##########################
# TASK 3 CONFIGURATION
# dQ vs pH Analysis
##########################

DEFAULT_RUN_Q_VS_PH = True  # Enable/disable Task 3 execution


##########################
# TASK 4 CONFIGURATION
# Onset Potential Analysis
# Rmenber set the ZPE in
#DEFAULT_E_ZPE_ST_FILE = [0, 0.334185, -0.00551, 0.308282]  # slab->ooh->o->oh
##########################

DEFAULT_RUN_ONSET_ANALYSIS = True     # Enable/disable Task 4 execution
DEFAULT_VOLTAGE_SCALE = "RHE"         # Voltage scale: "RHE" or "SHE" (must be RHE for onset calc)
DEFAULT_ONSET_PH_VALUES = [1, 13]     # pH values for onset potential calculations
DEFAULT_INCLUDE_O2 = False            # Include O2 in reaction intermediates

############################
# TASK 5 CONFIGURATION
# Single pH Point Analysis
# Rmenber set the ZPE in
#DEFAULT_E_ZPE_ST_FILE = [0, 0.334185, -0.00551, 0.308282]  # For file-based data,slab->ooh->o->oh
############################

DEFAULT_RUN_SINGLE_PH = True          # Enable/disable Task 5 execution
DEFAULT_ANALYSIS_PH_VALUES = [1,13]      # pH values for single point analysis

# Experimental RHE voltage ranges for OER/ORR
EXPERIMENTAL_RHE_RANGES = {
    "ORR": (-3, 2),     # Voltage range for Oxygen Reduction Reaction
    "OER": (1.2, 2.6)   # Voltage range for Oxygen Evolution Reaction
}


############################
# TASK 6 CONFIGURATION
# Heatmap Analysis
# Rmenber set the ZPE in
#DEFAULT_E_ZPE_ST_FILE = [0, 0.334185, -0.00551, 0.308282]  # For file-based data,slab->ooh->o->oh
############################

DEFAULT_RUN_HEATMAPS = True                 # Enable/disable Task 6 execution
DEFAULT_HEATMAP_PH_RANGE = np.linspace(0, 14, 150)     # pH range for heatmaps
DEFAULT_HEATMAP_VOLTAGE_RANGE = np.linspace(0, 1.5, 150)  # Voltage range for heatmaps


############################
# TASK 7 CONFIGURATION
# Reaction Pathway Diagrams
# Rmenber set the ZPE in
#DEFAULT_E_ZPE_ST_FILE = [0, 0.334185, -0.00551, 0.308282]  # For file-based data,slab->ooh->o->oh
############################

DEFAULT_RUN_STEP_DIAGRAMS = True  # Enable/disable Task 7 execution

# Conditions for pathway diagrams
DEFAULT_ORR_STEP_CONDITIONS = [
    {'label': 'pH=1, U=0 RHE (ORR)', 'pH': 1, 'U': 0},
    {'label': 'pH=13, U=0 RHE (ORR)', 'pH': 7, 'U': 0},
    #{'label': 'pH=1, U=0.4148 RHE (ORR)', 'pH': 1, 'U': 0.4148},#file
    {'label': 'pH=1, U=0.5077 RHE (ORR)', 'pH': 1, 'U': 0.5077}, # U= Onset potential
    #{'label': 'pH=13, U=0.6492 RHE (ORR)', 'pH': 13, 'U': 0.6492},#file
    {'label': 'pH=13, U=0.3344 RHE (ORR)', 'pH': 13, 'U': 0.3344} # U= Onset potential

]

DEFAULT_OER_STEP_CONDITIONS = [
    {'label': 'pH=1, U=0 RHE (OER)', 'pH': 1, 'U': 0},
    {'label': 'pH=13, U=0 RHE (OER)', 'pH': 13, 'U': 0},
   # {'label': 'pH=1, U=1.7016 RHE (OER)', 'pH': 1, 'U': 1.7016},#file
    #{'label': 'pH=13, U=1.6260 RHE (OER)', 'pH': 13, 'U': 1.6260},#file
    {'label': 'pH=1, U=1.7330 RHE (OER)', 'pH': 1, 'U': 1.7330}, #param
    {'label': 'pH=13, U=1.8771 RHE (OER)', 'pH': 13, 'U': 1.8771} #param
]


############################
# TASK 8 CONFIGURATION
# Microkinetics Parameters
# Rmenber set the ZPE in
#DEFAULT_E_ZPE_ST_FILE = [0, 0.334185, -0.00551, 0.308282]  # For file-based data,slab->ooh->o->oh
############################

DEFAULT_RUN_MICROKINETICS = True           # Enable/disable Task 8 execution
DEFAULT_MICROKINETICS_PH = [1, 13]         # pH values for microkinetics export
DEFAULT_MICROKINETICS_VOLTAGE = [0]        # Voltage values for microkinetics export


##############################
# TASK 9 CONFIGURATION
# Charge vs Potential Analysis
##############################

DEFAULT_RUN_QVSPOTENTIAL = True            # Enable/disable Task 9 execution
DEFAULT_DEFAULT_POTENTIAL_POINTS = [-1.23, 0.0, 1.23]  # Potential points for charge analysis
DEFAULT_CHARGE_PH = 1                      # Default pH for charge analysis


############################
# TASK 10 CONFIGURATION
# Adsorption Systems Analysis
############################

DEFAULT_RUN_OTHER_ADS_DATA = True          # Enable/disable Task 10 execution
ADS_SYSTEMS = ['*CO']                      # Adsorption systems to analyze
ADS_FILE_PATHS = ['CO.xlsx']               # Input files for adsorption systems
DEFAULT_U_MIN = -2                         # Min potential for adsorption plots
DEFAULT_U_MAX = 1                          # Max potential for adsorption plots


############################
# TASK 11-12 CONFIGURATION
# VASP Directory Creation and Data Extraction
############################

DEFAULT_RUN_VASP_DIR_CREATION = True       # Enable/disable Task 11 execution
DEFAULT_SUBMIT_VASP_JOBS = False           # Auto-submit VASP jobs after creation
DEFAULT_RUN_VASP_DATA_EXTRACTION = True    # Enable/disable Task 12 execution
DEFAULT_VASP_INPUT_FILES = ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'jy-sub_54_.vasp']  # Files to copy to VASP dirs
DEFAULT_VASP_SUBMIT_SCRIPT = 'jy-sub_54_.vasp'  # Submission script for VASP jobs


############################
# TASK 13 CONFIGURATION
# ORR/OER Microkinetic Modeling
############################

DEFAULT_RUN_ORR_OER_MICROKINETICS = True  # Enable/disable Task 13 execution

# Constants and parameters for microkinetic modeling
MICROKINETICS_CONFIG = {
    "acc": 50,                 # Numerical precision for sympy (int)
    "x_H2O": 1.0,              # Water mole fraction
    "x_O2": 2.34e-5,           # O2 mole fraction
    "T": 298.0,                # Temperature (K)
    "beta_i": [0, 0, 0.5, 0.5, 0.5, 0.5],  # Electron transfer coefficients
    "e": 1.602176634e-19,      # Elementary charge (C)
    "e_rho": 241 / 3 / 1000,   # Electron density factor (A/cm²/eV)
    "kBT": 8.617333e-5,        # Boltzmann constant * T (eV/K, pre-computed)
    "delta_Ga0i": [0.0, 0.0, 0.260, 0.260, 0.260, 0.260],  # Activation free energies (eV)
    "delta_Gi": [0.0, 0.0, -0.618389, -1.552131, -1.701575, -1.047905],  # Reaction free energies (eV)
    "Ai": [8.0e5, 1.0e8, 1.0e9, 1.0e9, 1.0e9, 1.0e9],  # Pre-exponential factors (s⁻¹)
    "current_max": 15,         # Max current for data truncation in plots (mA/cm²)
    "U_range_start": 0.0001,   # Start of U range (V vs RHE)
    "U_range_stop": 2.0,       # End of U range (V vs RHE)
    "U_range_step": 0.01       # Step size for U range (V)
}


############################
# TASK 14 CONFIGURATION
# OER/ORR Free Energy Diagrams
############################

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


#############################
# TASK 15 CONFIGURATION
# OER/ORR Activity Volcano Plot
#############################

DEFAULT_RUN_OER_ORR_VOLCANO = True   # Enable/disable Task 15

# Metal adsorption energies for volcano plot
# Format: {metal: (ΔG_OH, ΔG_O, ΔG_OOH)}
METAL_ADSORPTION_ENERGIES = {
    "Ti": (-0.93, -0.74, 2.5),
    "V": (-0.96, -0.99, 2.31),
    "Cr": (-0.48, -0.17, 2.89),
    "Mn": (-0.42, 0.47, 3.02),
    "Fe": (0.01, 0.75, 3.18),
    "Co": (0.23, 1.64, 3.39),
    "Ni": (0.13, 2.14, 3.39),
    "Cu": (0.75, 3.13, 4.07),
    "Zr": (-0.90, -0.39, 2.67),
    "Ru": (-0.30, 0.41, 2.83),
    "Rh": (0.78, 2.37, 3.9),
    "Pd": (0.41, 2.63, 3.65),
    "Ag": (0.97, 3.42, 4.25),
    "Ir": (0.58, 1.91, 3.71),
    "Pt": (0.19, 2.33, 3.56),
    "Au": (0.47, 2.55, 3.89)
}

# Plot configuration
VOLCANO_PLOT_RANGE = (-1.5, 3.0, -0.5, 3.5)  # (x_min, x_max, y_min, y_max)
OVERPTENTIAL_COLOR_RANGE = (0.32, 3.0)        # (vmin, vmax)

###########################
# EXECUTION FLAGS
# Master switches for tasks
###########################

# Enable/disable flags for all tasks
DEFAULT_RUN_Q_VS_PH = True          # Task 3
DEFAULT_RUN_ONSET_ANALYSIS = True   # Task 4
DEFAULT_RUN_SINGLE_PH = True        # Task 5
DEFAULT_RUN_HEATMAPS = True         # Task 6
DEFAULT_RUN_MICROKINETICS = True    # Task 8
DEFAULT_RUN_STEP_DIAGRAMS = True    # Task 7
DEFAULT_RUN_QVSPOTENTIAL = True     # Task 9
DEFAULT_RUN_OTHER_ADS_DATA = True   # Task 10
DEFAULT_RUN_VASP_DIR_CREATION = True     # Task 11
DEFAULT_RUN_VASP_DATA_EXTRACTION = True  # Task 12
DEFAULT_RUN_ORR_OER_MICROKINETICS = True  # Task 13
