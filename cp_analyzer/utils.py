# utils.py
import os
import re
from unidecode import unidecode
import numpy as np
import pandas as pd

def create_output_directory(path):
    """Create output directory and return path"""
    os.makedirs(path, exist_ok=True)
    return path

def safe_filename(name):
    """Create filesystem-safe name (remove or translate special characters)"""
    safe = re.sub(r'[^\w_. -]', '', name)
    if not safe or len(safe) < 1:
        safe = unidecode(name)
        safe = re.sub(r'[^\w_. -]', '', safe)
    if not safe or len(safe) < 1:
        safe = "system_" + str(hash(name))
    return safe

def U_rhe_to_U_she(U_vs_rhe, pH):
    """Convert RHE voltage to SHE voltage"""
    return U_vs_rhe - 0.0592 * pH

def U_she_to_U_rhe(U_vs_she, pH):
    """Convert SHE voltage to RHE voltage"""
    return U_vs_she + 0.0592 * pH

def get_experimental_voltage_range(pH, reaction_type, voltage_scale, rhe_ranges):
    rhe_range = rhe_ranges[reaction_type]
    if voltage_scale == "RHE":
        return rhe_range
    else:
        min_voltage = rhe_range[0] - 0.0592 * pH
        max_voltage = rhe_range[1] - 0.0592 * pH
        return (min_voltage, max_voltage)

def E_free_vs_SHE(E_fermis, E_Vacums, E_potentials, E_DFTs, Charges):
    E_fermi_vs_Vacums = E_fermis - E_Vacums
    E_potential_vs_Vacums = E_potentials - E_Vacums
    E_applied_potentials = E_potential_vs_Vacums
    E_potential_0s = np.array([E_applied_potentials[0]] * len(E_applied_potentials))
    E_corrs = (E_applied_potentials + E_potential_0s) * Charges / 2
    E_q_Wf = E_fermi_vs_Vacums * Charges
    E_free_U = E_DFTs + E_corrs - E_q_Wf
    E_U_vs_SHE = -4.6 - E_fermi_vs_Vacums
    return E_free_U, E_U_vs_SHE

def quadratic_function_vs_she(E_U_vs_SHE, C, U_pzc, E0):
    return -0.5 * C * (E_U_vs_SHE - U_pzc) ** 2 + E0

def read_excel_columns_by_index(file_path):
    try:
        df = pd.read_excel(file_path)
        return [df.iloc[:, i].tolist() for i in range(len(df.columns))]
    except Exception as e:
        print(f"Error reading file: {e}")
        return None