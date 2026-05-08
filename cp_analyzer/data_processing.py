# data_processing.py
import numpy as np
from scipy.optimize import curve_fit

# 包内导入
from cp_analyzer.utils import quadratic_function_vs_she, read_excel_columns_by_index, E_free_vs_SHE

import os
import re
from unidecode import unidecode
import pandas as pd
def get_fit_parameters_from_data(data):
    """Get fitting parameters from data structure, including R²"""
    C_list, U_pzc_list, E0_list, R2_list = [], [], [], []
    systems = list(data.keys())
    for system in systems:
        U = np.array(data[system]["U"])
        E = np.array(data[system]["E"])
        if len(U) < 3 or len(E) < 3:  # Ensure sufficient data points
            print(f"Warning: System {system} has insufficient data points, skipping fit")
            continue

        try:
            # Perform fitting
            params, cov = curve_fit(quadratic_function_vs_she, U, E, p0=[1.0, 0.0, 0.0])
            C, U_pzc, E0 = params

            # Calculate R²
            E_pred = quadratic_function_vs_she(U, C, U_pzc, E0)
            ss_res = np.sum((E - E_pred) ** 2)
            ss_tot = np.sum((E - np.mean(E)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else float('nan')

            C_list.append(C)
            U_pzc_list.append(U_pzc)
            E0_list.append(E0)
            R2_list.append(r_squared)

            print(f"{system}: C={C:.4f}, U_pzc={U_pzc:.4f}, E0={E0:.4f}, R²={r_squared:.4f}")
        except Exception as e:
            print(f"Error fitting {system}: {e}")
            continue
    return C_list, U_pzc_list, E0_list, R2_list
def fit_ads_quadratic(E_U_vs_SHE, E_free_U):
    """Fit quadratic function for adsorption energy vs SHE"""

    def quadratic_func(U, C, U_pzc, E0):
        return -0.5 * C * (U - U_pzc) ** 2 + E0

    params, cov = curve_fit(quadratic_func, E_U_vs_SHE, E_free_U, p0=[1.0, 0.0, 0.0])
    C, U_pzc, E0 = params

    # Calculate R²
    E_pred = quadratic_func(E_U_vs_SHE, *params)
    ss_res = np.sum((E_free_U - E_pred) ** 2)
    ss_tot = np.sum((E_free_U - np.mean(E_free_U)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else float('nan')

    return C, U_pzc, E0, r_squared

def get_fit_parameters_from_files(file_paths):
    """
    Get fitting parameters (C, U_pzc, E0) from Excel files, including R².
    Args:
        file_paths (list): List of paths to Excel files containing electrochemical data.
    Returns:
        tuple: Lists of C, U_pzc, E0, and R² values for each file.
    """
    C_list, U_pzc_list, E0_list, R2_list = [], [], [], []

    for file_path in file_paths:
        # Read data from Excel
        data_columns = read_excel_columns_by_index(file_path)
        if not data_columns or len(data_columns) < 5:
            print(f"File {file_path} has insufficient columns, skipping")
            continue

        # Extract data
        E_fermis = np.array(data_columns[0])
        E_Vacums = np.array(data_columns[1])
        E_potentials = np.array(data_columns[2])
        E_DFTs = np.array(data_columns[3])
        Charges = np.array(data_columns[4])

        # Validate input data
        if not all(len(arr) > 0 for arr in [E_fermis, E_Vacums, E_potentials, E_DFTs, Charges]):
            print(f"File {file_path} has empty data arrays, skipping")
            continue

        # Calculate free energy
        try:
            E_free_U, E_U_vs_SHE = E_free_vs_SHE(
                E_fermis, E_Vacums, E_potentials, E_DFTs, Charges
            )

            # Remove invalid points
            mask = np.isfinite(E_U_vs_SHE) & np.isfinite(E_free_U)
            E_U_vs_SHE = E_U_vs_SHE[mask]
            E_free_U = E_free_U[mask]

            if len(E_U_vs_SHE) < 3:
                print(f"File {file_path} has insufficient valid data points ({len(E_U_vs_SHE)}), skipping")
                continue

            # Print data for debugging
            print(f"File {file_path}:")
            print(f"  E_U_vs_SHE: {E_U_vs_SHE}")
            print(f"  E_free_U: {E_free_U}")

        except Exception as e:
            print(f"Error calculating free energy for {file_path}: {e}")
            continue

        try:
            # Dynamic initial guess
            C_init = 1.0  # Typical capacitance-like coefficient
            U_pzc_init = np.mean(E_U_vs_SHE) if len(E_U_vs_SHE) > 0 else 0.0  # Center of potential range
            E0_init = np.mean(E_free_U) if len(E_free_U) > 0 else 0.0  # Center of energy range

            # Parameter bounds: C > 0, U_pzc within reasonable voltage, E0 within reasonable energy
            bounds = ([0, -5, -1000], [10, 5, 1000])  # Adjust based on expected ranges

            # Perform fitting with bounds and increased maxfev
            params, cov = curve_fit(
                quadratic_function_vs_she,
                E_U_vs_SHE,
                E_free_U,
                p0=[C_init, U_pzc_init, E0_init],
                bounds=bounds,
                maxfev=10000,  # Increase max function evaluations
                method='trf'  # Trust Region Reflective for robustness
            )
            C, U_pzc, E0 = params

            # Calculate R²
            E_pred = quadratic_function_vs_she(E_U_vs_SHE, C, U_pzc, E0)
            ss_res = np.sum((E_free_U - E_pred) ** 2)
            ss_tot = np.sum((E_free_U - np.mean(E_free_U)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else float('nan')

            # Store parameters
            C_list.append(C)
            U_pzc_list.append(U_pzc)
            E0_list.append(E0)
            R2_list.append(r_squared)

            file_name = os.path.splitext(os.path.basename(file_path))[0]
            print(f"{file_name}: C={C:.4f}, U_pzc={U_pzc:.4f}, E0={E0:.4f}, R²={r_squared:.4f}")

        except Exception as e:
            print(f"Error fitting {file_path}: {e}")
            continue

    return C_list, U_pzc_list, E0_list, R2_list
def create_output_directory(path):
    """Create output directory and return path"""
    os.makedirs(path, exist_ok=True)
    return path

def get_fit_params_from_literature_params(literature_params):
    """Get fitting parameters directly from literature parameters"""
    C_list, U_pzc_list, E0_list = [], [], []

    # Fixed species order: must be slab, *OOH, *O, *OH
    species_order = ["slab", "*OOH", "*O", "*OH"]

    for spec in species_order:
        if spec in literature_params:
            C_list.append(literature_params[spec]["C"])
            U_pzc_list.append(literature_params[spec]["U_pzc"])
            E0_list.append(literature_params[spec]["E0"])
        else:
            print(f"Warning: {spec} not found in literature parameters, setting default to 0")
            C_list.append(0)
            U_pzc_list.append(0)
            E0_list.append(0)
    return C_list, U_pzc_list, E0_list
def calculate_ads_energy(U, C, U_pzc, E0, E_ref=0.0):
    """Calculate adsorption energy at given U (SHE)"""
    E_ads = -0.5 * C * (U - U_pzc) ** 2 + E0 - E_ref  # ΔE = E(U) - E_ref, assuming slab is reference
    return E_ads
def safe_filename(name):
    """Create filesystem-safe name (remove or translate special characters)"""
    # 第一步：移除非法字符（保留字母数字、下划线、点、空格和短横线）
    safe = re.sub(r'[^\w_. -]', '', name)

    # 第二步：如果结果为空或长度为0，使用unidecode转换
    if not safe or len(safe) < 1:  # 修改的关键：允许长度≥1
        safe = unidecode(name)  # 转换Unicode到ASCII
        safe = re.sub(r'[^\w_. -]', '', safe)  # 再次移除非法字符

    # 第三步：确保最低长度（至少1字符）
    if not safe or len(safe) < 1:
        safe = "system_" + str(hash(name))
    return safe

def save_ads_data(output_dir, ads_systems, fit_params_list, U_range):
    """Save adsorption energy data and fitting parameters"""
    data_dir = create_output_directory(os.path.join(output_dir, "ads_data"))
    for i, system in enumerate(ads_systems):
        C, U_pzc, E0, r2 = fit_params_list[i]
        E_ads = [calculate_ads_energy(u, C, U_pzc, E0) for u in U_range]
        df = pd.DataFrame({"U": U_range, "Delta_E": E_ads})
        safe_system = safe_filename(system)  # Use safe_filename to handle invalid characters
        df.to_csv(os.path.join(data_dir, f"{safe_system}_ads_energy.csv"), index=False)

        # Save fitting parameters
        param_file = os.path.join(data_dir, f"{safe_system}_fit_params.txt")
        with open(param_file, 'w') as f:
            f.write(f"System: {system}\n")
            f.write(f"C: {C:.6f}\n")
            f.write(f"U_pzc: {U_pzc:.6f}\n")
            f.write(f"E0: {E0:.6f}\n")
            f.write(f"R2: {r2:.6f}\n")

        print(f"Fitting parameters for {system}: C={C:.6f}, U_pzc={U_pzc:.6f}, E0={E0:.6f}, R2={r2:.6f}")
def read_excel_columns_by_index(file_path):
    """Read data columns from Excel file"""
    try:
        df = pd.read_excel(file_path)
        return [df.iloc[:, i].tolist() for i in range(len(df.columns))]
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def Q_tran_vs_pH(file_paths=None, Urhe=0):
    """Calculate Q vs pH relation"""
    Q_4 = []
    pH_4 = []
    for file_path in file_paths:
        data_columns = read_excel_columns_by_index(file_path)
        if not data_columns or len(data_columns) < 6:
            continue
        E_free_U, E_U_vs_SHE = E_free_vs_SHE(
            np.array(data_columns[0]), np.array(data_columns[1]),
            np.array(data_columns[2]), np.array(data_columns[3]),
            np.array(data_columns[4]))

        # Calculate pH value
        pH = (Urhe - E_U_vs_SHE) / 0.0592
        Q = np.array(data_columns[5])  # 6th column is charge
        Q_4.append(Q)
        pH_4.append(pH)
    return Q_4, pH_4