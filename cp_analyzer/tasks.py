# tasks.py
from __future__ import annotations

# 包内正确导入
from cp_analyzer.config_default import *
from cp_analyzer.data_processing import *
from cp_analyzer.plotting import *
from cp_analyzer.calculations import *
from cp_analyzer.utils import *

import shutil
import subprocess
import pandas as pd
from scipy.integrate import simpson
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from typing import List, Tuple
def task1(args):
    print("=" * 80)
    print("TASK 1: Obtaining Fitting Parameters".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_1_fitting_parameters"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    if args.data_source == "literature_params":
        C_list, U_pzc_list, E0_list = get_fit_params_from_literature_params(LITERATURE_PARAMS)
        R2_list = [None] * len(C_list)
        data_to_use = LITERATURE_PARAMS
    elif args.data_source == "literature":
        C_list, U_pzc_list, E0_list, R2_list = get_fit_parameters_from_data(LITERATURE_DATA)
        data_to_use = LITERATURE_DATA
    else:
        file_paths = args.file_paths
        if not file_paths:
            print("ERROR: No file paths provided")
            return None, None, None, None
        C_list, U_pzc_list, E0_list, R2_list = get_fit_parameters_from_files(file_paths)
        data_to_use = None

    if C_list and U_pzc_list and E0_list:
        if args.include_o2 and len(C_list) >= 5:
            names = ["Slab", "O2", "*OOH", "*O", "*OH"]
        elif len(C_list) >= 4:
            names = ["Slab", "*OOH", "*O", "*OH"]
        else:
            names = [f"System {i + 1}" for i in range(len(C_list))]

        print("\nSummary of fitting results:")
        print(f"{'System':<16} | {'Coefficient (C)':>15} | {'U_pzc (V)':>10} | {'E0 (eV)':>10} | {'R^2':>8}")
        print("-" * 60)

        for i in range(min(len(C_list), len(names))):
            r2_value = R2_list[i] if i < len(R2_list) else None
            if isinstance(r2_value, (int, float)) and not np.isnan(r2_value):
                r2_str = f"{r2_value:>8.4f}"
            else:
                r2_str = "N/A".rjust(8)

            print(f"{names[i]:<16} | {C_list[i]:15.6f} | {U_pzc_list[i]:10.6f} | "
                  f"{E0_list[i]:10.6f} | {r2_str}")
    # ✅ Optional: save fitting parameters to file (Task 1 output)
    if C_list and U_pzc_list and E0_list:
        save_path = os.path.join(task_output_dir, "fitting_parameters.csv")
        with open(save_path, "w", encoding="utf-8") as f:
            f.write("System,C,U_pzc,E0,R2\n")
            for i in range(len(C_list)):
                system_name = names[i] if i < len(names) else f"System_{i+1}"
                r2_value = R2_list[i] if i < len(R2_list) else None
                r2_str = f"{r2_value:.6f}" if isinstance(r2_value, (int, float)) else "N/A"
                f.write(
                    f"{system_name},{C_list[i]:.8f},{U_pzc_list[i]:.8f},"
                    f"{E0_list[i]:.8f},{r2_str}\n"
                )

        print(f"\n✓ Fitting parameters saved to: {save_path}")

    return C_list, U_pzc_list, E0_list, R2_list, data_to_use

def task2(args, C_list, U_pzc_list, E0_list, R2_list, data_to_use):
    task2_output_dir = os.path.join(
        args.output_dir,
        "TASK_2_free_energy_fits"
    )
    os.makedirs(task2_output_dir, exist_ok=True)

    print("=" * 80)
    print("TASK 2: Plotting Free Energy vs SHE Potential".center(80))
    print("=" * 80)

    file_paths = args.file_paths if args.data_source == "file" else []
    plot_E_free_vs_U_she(task2_output_dir, args.data_source, data_to_use, file_paths, R2_list)

def task3(args, file_paths):
    """
    TASK 3: dQ vs pH Analysis
    Analyze charge change (dQ) of intermediates as a function of pH.
    """

    import os
    import pandas as pd

    # ============================================================
    # Task header (always print)
    # ============================================================
    print("=" * 80)
    print("TASK 3: dQ vs pH Analysis".center(80))
    print("=" * 80)

    # ============================================================
    # Check whether task is enabled
    # ============================================================
    if not args.run_q_vs_ph:
        print(
            "TASK 3 skipped: --run-q-vs-ph not specified.\n"
            "Please run:\n"
            "  python main.py --task 3 --run-q-vs-ph "
            "--data-source file --file-paths your_data.xlsx"
        )
        return False

    # ============================================================
    # Check data source
    # ============================================================
    if args.data_source != "file":
        print(
            f"TASK 3 skipped: --data-source must be 'file'.\n"
            f"Current value: {args.data_source}"
        )
        return False

    # ============================================================
    # Check input files
    # ============================================================
    if not file_paths:
        print(
            "TASK 3 skipped: no input files provided.\n"
            "Please specify:\n"
            "  --file-paths your_data.xlsx"
        )
        return False

    # ============================================================
    # Load data and explicitly search for dQ column
    # ============================================================
    valid_dq_data = False

    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"Warning: file not found, skipping: {file_path}")
            continue

        try:
            df = pd.read_excel(file_path)

            # ----------------------------------------------------
            # Explicitly require a dQ column (strict check)
            # ----------------------------------------------------
            dq_columns = [
                col for col in df.columns
                if str(col).strip().lower() in {
                    "dq",
                    "delta_q",
                    "charge_change"
                }
            ]

            if not dq_columns:
                # This file does not meet the requirement
                continue

            dq_column = dq_columns[0]

            # Check whether dQ column has valid numeric values
            dq_values = pd.to_numeric(df[dq_column], errors="coerce").dropna()

            if not dq_values.empty:
                valid_dq_data = True
                # ---- your original dQ vs pH processing logic here ----
                # e.g. collect data, plot, save results

        except Exception as e:
            print(f"Warning: failed to read {file_path}: {e}")

    # ============================================================
    # No valid dQ data → print ORIGINAL warning (guaranteed)
    # ============================================================
    if not valid_dq_data:
        print("=" * 80)
        print(
            "Warning: No valid dQ/pH data, skipping dQ vs pH analysis,\n"
            "you need input dQ(the charge change)\n"
            "between intermediates and adsorption center in *.xlsx at the last column!!"
        )
        print("=" * 80)
        return False

    # ============================================================
    # If reached here, task executed successfully
    # ============================================================
    print("TASK 3 completed: dQ vs pH analysis finished successfully.")
    return True


def task4(args, C_list, U_pzc_list, E0_list, E_zpe_st):
    if not args.run_onset_analysis:
        return
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_4_onset_potential"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    print("=" * 80)
    print("TASK 4: Onset Potential Analysis".center(80))
    print("=" * 80)

    if args.voltage_scale != "RHE":
        print("ERROR: Onset potential calculation requires RHE voltage scale")
        return

    calculate_and_output_onset_potentials(
        task_output_dir,
        args.onset_ph_values,
        C_list, U_pzc_list, E0_list, E_zpe_st,
        args.include_o2,
        args.voltage_scale
    )


def task5(args, C_list, U_pzc_list, E0_list, E_zpe_st):
    if not args.run_single_ph:
        return

    print("=" * 80)
    print("TASK 5: Single pH Point Analysis".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_5_single_pH_analysis"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    pH_dir = create_output_directory(os.path.join(task_output_dir, "pH_analysis"))

    for pH in args.analysis_ph_values:
        current_dir = create_output_directory(os.path.join(pH_dir, f"pH_{pH}"))

        for reaction_type in ["ORR", "OER"]:
            rhe_min, rhe_max = EXPERIMENTAL_RHE_RANGES[reaction_type]

            if args.voltage_scale == "RHE":
                voltage_range = np.linspace(rhe_min, rhe_max, 50)
            else:
                she_min = U_rhe_to_U_she(rhe_min, pH)
                she_max = U_rhe_to_U_she(rhe_max, pH)
                voltage_range = np.linspace(she_min, she_max, 50)

            results = calculate_for_voltage_range(voltage_range, pH, C_list, U_pzc_list, E0_list, E_zpe_st, reaction_type, args.include_o2, args.voltage_scale)

            plot_dG_diagrams(current_dir, results, pH, reaction_type, args.voltage_scale)
            save_single_pH_results(current_dir, results, pH, reaction_type, args.voltage_scale)

def task6(args, C_list, U_pzc_list, E0_list, E_zpe_st):
    if not args.run_heatmaps:
        return
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_6_heatmap_analysis"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    print("=" * 80)
    print("TASK 6: Heatmap Analysis".center(80))
    print("=" * 80)

    if args.voltage_scale != "RHE":
        print("ERROR: Heatmap generation requires RHE voltage scale")
        return

    for reaction_type in ["ORR", "OER"]:
        grid_data = calculate_dG_grid(args.heatmap_ph_range, args.heatmap_voltage_range, C_list, U_pzc_list, E0_list, E_zpe_st, reaction_type, args.include_o2, args.voltage_scale)
        save_dG_grid(grid_data, task_output_dir)
        plot_dG_heatmaps(grid_data, task_output_dir)


def task7(args, C_list, U_pzc_list, E0_list, E_zpe_st):
    if not args.run_step_diagrams:
        return

    print("=" * 80)
    print("TASK 7: Reaction Pathway Diagram Analysis".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_7_reaction_pathway"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    if args.voltage_scale != "RHE":
        print("ERROR: Step diagram generation requires RHE voltage scale")
        return

    plot_reaction_pathway(task_output_dir, DEFAULT_ORR_STEP_CONDITIONS, C_list, U_pzc_list, E0_list, E_zpe_st, "ORR", args.include_o2, args.voltage_scale)
    plot_reaction_pathway(task_output_dir, DEFAULT_OER_STEP_CONDITIONS, C_list, U_pzc_list, E0_list, E_zpe_st, "OER", args.include_o2, args.voltage_scale)

def task8(args, C_list, U_pzc_list, E0_list, E_zpe_st):
    if not args.run_microkinetics:
        return

    print("=" * 80)
    print("TASK 8: Microkinetics Parameter Export".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_8_microkinetics_parameters"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    if args.voltage_scale != "RHE":
        print("ERROR: Microkinetics parameter generation requires RHE voltage scale")
        return

    save_microkinetics_parameters(task_output_dir, DEFAULT_MICROKINETICS_PH, DEFAULT_MICROKINETICS_VOLTAGE, C_list, U_pzc_list, E0_list, E_zpe_st, args.include_o2, args.voltage_scale)

def task9(args, file_paths):
    if not args.run_qvspotential:
        return

    print("=" * 80)
    print("TASK 9: Charge vs Potential Analysis".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_9_charge_vs_potential"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    if args.data_source == "file" and file_paths:
        plot_charge_vs_potential(task_output_dir, file_paths, plot_type="combined", scale_type=args.voltage_scale, charge_pH=args.charge_ph, DEFAULT_POTENTIAL_POINTS=args.default_potential_points)

def task10(args, file_paths):
    if not args.run_other_ads_data:
        return

    print("=" * 80)
    print("TASK 10: Generic Adsorption Systems Analysis".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_10_adsorption_analysis"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    # Use configuration from config.py
    fit_params_list = []
    names = ADS_SYSTEMS  # System names from config

    for i, path in enumerate(ADS_FILE_PATHS):  # Updated to use ADS_FILE_PATHS from config
        data_columns = read_excel_columns_by_index(path)
        if not data_columns or len(data_columns) < 5:
            print(f"Warning: Invalid or insufficient data in {path}. Skipping.")
            continue
        E_free_U, E_U_vs_SHE = E_free_vs_SHE(
            np.array(data_columns[0]), np.array(data_columns[1]),
            np.array(data_columns[2]), np.array(data_columns[3]),
            np.array(data_columns[4])
        )
        C, U_pzc, E0, r2 = fit_ads_quadratic(E_U_vs_SHE, E_free_U)
        fit_params_list.append((C, U_pzc, E0, r2))

    # Print summary of fitting results (like task1)
    if fit_params_list:
        print("\nSummary of fitting results:")
        print(f"{'System':<16} | {'Coefficient (C)':>15} | {'U_pzc (V)':>10} | {'E0 (eV)':>10} | {'R^2':>8}")
        print("-" * 60)
        for i in range(min(len(fit_params_list), len(names))):
            C, U_pzc, E0, r2 = fit_params_list[i]
            if isinstance(r2, (int, float)) and not np.isnan(r2):
                r2_str = f"{r2:>8.4f}"
            else:
                r2_str = "N/A".rjust(8)
            print(f"{names[i]:<16} | {C:15.6f} | {U_pzc:10.6f} | "
                  f"{E0:10.6f} | {r2_str}")

    U_range = np.linspace(args.u_min, args.u_max, 100)

    # Plot and save data (like task2: plot + save points)
    plot_ads_energy_vs_U(task_output_dir, ADS_SYSTEMS, fit_params_list, "SHE", pH=0, U_range=U_range)  # Updated to use ADS_SYSTEMS
    save_ads_data(task_output_dir, ADS_SYSTEMS, fit_params_list, U_range)  # Saves fitting params + computed energies to CSV

    plot_ads_energy_vs_U(task_output_dir, ADS_SYSTEMS, fit_params_list, "RHE", pH=0, U_range=U_range)  # Updated to use ADS_SYSTEMS


def task11(args, vasp_input_files=None):
    """
    Task 11: Create VASP calculation directories with modified NELECT values and optionally submit jobs.

    This task creates directories for VASP calculations with different NELECT increments relative to the neutral system.
    It copies VASP input files (INCAR, KPOINTS, POSCAR, POTCAR, and submit script) to each directory and modifies
    the NELECT value in INCAR based on the specified increments. Optionally submits jobs using sbatch.
    """
    if not args.run_vasp_dir_creation:
        return

    print("=" * 80)
    print("TASK 11: CREATE VASP CALCULATION DIRECTORIES".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_11_vasp_directories"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    # Define default VASP input files if not provided
    if vasp_input_files is None:
        vasp_input_files = ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', args.vasp_submit_script]

    # Validate input files
    missing_files = [f for f in vasp_input_files if not os.path.exists(f)]
    if missing_files:
        print(f"ERROR: The following required VASP input files are missing: {missing_files}")
        print("Task 11 aborted. Please ensure all VASP input files are in the current directory.")
        return

    # Read NELECT from INCAR file
    try:
        nelect = None
        with open('INCAR', 'r', encoding='utf-8') as incar_file:
            for line in incar_file:
                line = line.strip()
                if line.startswith('NELECT =') or ('NELECT' in line and '=' in line):
                    try:
                        # Handle different formats: NELECT = 100.0 or NELECT=100.0
                        nelect_str = line.split('=')[1].strip()
                        nelect = float(nelect_str)
                        break
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse NELECT from line '{line}': {e}")
                        continue
        if nelect is None:
            print("ERROR: NELECT tag not found or could not be parsed from INCAR file")
            print("Please ensure INCAR contains a valid 'NELECT = X.X' line.")
            return
    except FileNotFoundError:
        print("ERROR: INCAR file not found in current directory")
        return
    except Exception as e:
        print(f"ERROR: Failed to read INCAR file: {e}")
        return

    print(f"Neutral system NELECT = {nelect}")

    # Define increments and corresponding folder names (matching original script)
    increments = [0, -0.5, -1, -1.5, -2, 0.5, 1, 1.5, 2]
    folder_names = ['e0', 'e-0.5', 'e-1', 'e-1.5', 'e-2', 'e+0.5', 'e+1', 'e+1.5', 'e+2']

    # Create output directory for VASP calculations
    vasp_base_dir = create_output_directory(
        os.path.join(task_output_dir, "vasp_calculations")
    )

    # Track successful creations and submissions
    created_dirs = 0
    submitted_jobs = 0

    # Process each increment
    for folder_name, increment in zip(folder_names, increments):
        folder_path = os.path.join(vasp_base_dir, folder_name)
        os.makedirs(folder_path, exist_ok=True)

        # Copy all specified input files to the new directory
        copied_files = []
        for file_name in vasp_input_files:
            try:
                shutil.copy(file_name, folder_path)
                copied_files.append(file_name)
            except Exception as e:
                print(f"Warning: Failed to copy {file_name} to {folder_name}: {e}")
                continue

        if len(copied_files) < len(vasp_input_files):
            print(f"Warning: Incomplete file copy for {folder_name}. Skipping INCAR modification.")
            continue

        # Modify INCAR file if increment is non-zero
        if increment != 0:
            new_nelect = nelect + increment
            incar_path = os.path.join(folder_path, 'INCAR')
            try:
                # Read existing INCAR lines
                with open(incar_path, 'r', encoding='utf-8') as incar_file:
                    lines = incar_file.readlines()

                # Write back with modified NELECT
                with open(incar_path, 'w', encoding='utf-8') as incar_file:
                    for line in lines:
                        if line.strip().startswith('NELECT =') or ('NELECT' in line and '=' in line):
                            # Replace the NELECT line
                            incar_file.write(f'NELECT = {new_nelect}\n')
                        else:
                            incar_file.write(line)

                print(f"  Created {folder_name}: NELECT = {new_nelect} (increment: {increment:+.1f})")
            except Exception as e:
                print(f"ERROR: Failed to modify INCAR in {folder_name}: {e}")
                continue

        else:
            print(f"  Created {folder_name}: NELECT = {nelect} (neutral)")

        created_dirs += 1

        # Optionally submit job if enabled
        if args.submit_vasp_jobs:
            submit_script_name = args.vasp_submit_script
            submit_script_path = os.path.join(folder_path, submit_script_name)
            if os.path.exists(submit_script_path):
                try:
                    # Change to the folder directory for submission
                    original_cwd = os.getcwd()
                    os.chdir(folder_path)

                    # Submit the job using sbatch
                    result = subprocess.run(['sbatch', submit_script_name],
                                            capture_output=True, text=True, check=True)
                    print(f"  Submitted job in {folder_name}: {result.stdout.strip()}")
                    submitted_jobs += 1

                    # Return to original directory
                    os.chdir(original_cwd)

                except subprocess.CalledProcessError as e:
                    print(f"  ERROR submitting job in {folder_name}: {e.stderr}")
                except FileNotFoundError:
                    print(f"  ERROR: sbatch command not found. Please run jobs manually.")
                except Exception as e:
                    print(f"  ERROR during job submission in {folder_name}: {e}")
                    if 'original_cwd' in locals():
                        os.chdir(original_cwd)
            else:
                print(f"  Warning: Submit script '{submit_script_name}' not found in {folder_name}")

    # Summary
    print(f"\n✓ Task 11 completed successfully!")
    print(f"  Created {created_dirs} directories in: {vasp_base_dir}")
    if args.submit_vasp_jobs:
        print(f"  Submitted {submitted_jobs} jobs automatically")
        if submitted_jobs < created_dirs:
            print(f"  Note: {created_dirs - submitted_jobs} jobs were not submitted (check warnings above)")
    else:
        print(f"  Note: Job submission disabled (--submit-vasp-jobs flag not used)")
        print(f"  To submit jobs manually, cd into each directory and run: sbatch {args.vasp_submit_script}")

    print(f"\nNELECT values relative to neutral system (NELECT = {nelect}):")
    print("  " + " | ".join([f"{folder}: {nelect + inc:+.1f}" for folder, inc in zip(folder_names, increments)]))
    print("=" * 80)



def task12(args):
    """
    Task 12: Extract VASP data from directories and save to Excel.

    This task reads VASP output files (LOCPOT, OUTCAR) from directories created
    for different charge states (e.g., e0, e-0.5, e-1, etc.), calculates
    relevant electrochemical properties (Fermi energy, vacuum levels, average
    potential, DFT energy), and saves the results to an Excel file.
    """
    if not args.run_vasp_data_extraction:
        return

    print("=" * 80)
    print("TASK 12: EXTRACT VASP DATA TO EXCEL".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_12_vasp_data_extraction"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    # Define charges and corresponding folder names (matching Task 11)
    charges = [0, -0.5, -1, -1.5, -2, 0.5, 1, 1.5, 2]
    folder_names = ['e0', 'e-0.5', 'e-1', 'e-1.5', 'e-2', 'e+0.5', 'e+1', 'e+1.5', 'e+2']
    columns = ['E_fermis', 'top_E_Vacums', 'E_potentials', 'E_DFTs', 'Charges']

    # Base directory for VASP calculations
    vasp_base_dir = os.path.join(
        args.output_dir,
        "TASK_11_vasp_directories",
        "vasp_calculations"
    )

    if not os.path.exists(vasp_base_dir):
        print(f"ERROR: VASP calculations directory not found at {vasp_base_dir}")
        print("Task 12 aborted. Please ensure Task 11 has been run to create directories.")
        return

    data_list = []

    # Process each folder
    for folder_name, charge in zip(folder_names, charges):
        folder_path = os.path.join(vasp_base_dir, folder_name)
        print(f"Extracting data from: {folder_name}...")

        if not os.path.exists(folder_path):
            print(f"  Warning: Directory {folder_name} not found. Skipping.")
            continue

        try:
            os.chdir(folder_path)

            # Helper function: Read LOCPOT for z-length and plane averages
            def read_z_locpot():
                try:
                    with open('LOCPOT', 'r', encoding='utf-8') as file:
                        content = file.readlines()

                    # Extract lattice vector for c-direction
                    coords_c = list(map(float, content[4].strip().split()))
                    length_c = np.linalg.norm(coords_c)

                    # Find grid information
                    data_start_index = int(sum(map(float, content[6].split()))) + 10
                    grid_info = list(map(int, content[data_start_index - 1].strip().split()))
                    total_num = int((grid_info[0] * grid_info[1] * grid_info[2]) / 5)

                    # Extract potential data
                    data = []
                    for line in content[data_start_index:data_start_index + total_num]:
                        values = list(map(float, line.strip().split()))
                        if len(values) != 5:
                            print(f"  Warning: Invalid data line in LOCPOT at {folder_name}")
                            continue
                        data.extend(values)

                    # Reshape to 3D array
                    data = np.array(data).reshape((grid_info[2], grid_info[0], grid_info[1]))
                    z_lengths = np.linspace(0, length_c, grid_info[2])
                    plane_averages = np.mean(data, axis=(1, 2))

                    return length_c, z_lengths, plane_averages
                except Exception as e:
                    print(f"  ERROR: Failed to read LOCPOT in {folder_name}: {e}")
                    return None, None, None

            # Helper function: Calculate average potential
            def calculate_average_potential(x_values, y_values):
                if x_values is None or y_values is None:
                    return None
                try:
                    total_area = simpson(y=y_values, x=x_values)
                    average_potential = total_area / x_values[-1]
                    print(f"  E_potential = {average_potential:.4f}")
                    return average_potential
                except Exception as e:
                    print(f"  ERROR: Failed to calculate average potential in {folder_name}: {e}")
                    return None

            # Helper function: Calculate vacuum averages
            def calculate_vacuum_averages(y_values, threshold=0.01, window_size=5):
                if y_values is None:
                    return None, None
                try:
                    start_idx = 0
                    while start_idx < len(y_values) - window_size:
                        if max(y_values[start_idx:start_idx + window_size]) - min(
                                y_values[start_idx:start_idx + window_size]) <= threshold:
                            start_idx += 1
                        else:
                            break

                    end_idx = len(y_values) - 1
                    while end_idx > window_size:
                        if max(y_values[end_idx - window_size:end_idx]) - min(
                                y_values[end_idx - window_size:end_idx]) <= threshold:
                            end_idx -= 1
                        else:
                            break

                    front_vacuum_average = np.mean(y_values[:start_idx + 1]) if start_idx > 0 else 0
                    end_vacuum_average = np.mean(y_values[end_idx:]) if end_idx < len(y_values) - 1 else 0
                    print(f"  top_Vacuum_level = {front_vacuum_average:.4f}")
                    print(f"  bottom_Vacuum_level = {end_vacuum_average:.4f}")
                    return front_vacuum_average, end_vacuum_average
                except Exception as e:
                    print(f"  ERROR: Failed to calculate vacuum averages in {folder_name}: {e}")
                    return None, None

            # Helper function: Read Fermi energy from OUTCAR
            def read_fermi_energy():
                try:
                    with open('OUTCAR', 'r', encoding='utf-8') as file:
                        for line in file:
                            if 'E-fermi' in line:
                                e_fermi = float(line.split()[2])
                                print(f"  E-fermi = {e_fermi:.4f}")
                                return e_fermi
                    print(f"  Warning: E-fermi not found in OUTCAR for {folder_name}")
                    return None
                except Exception as e:
                    print(f"  ERROR: Failed to read OUTCAR in {folder_name}: {e}")
                    return None

            # Helper function: Extract TOTEN from OUTCAR
            def extract_toten_from_outcar():
                try:
                    toten_values = []
                    with open('OUTCAR', 'r', encoding='utf-8') as file:
                        for line in file:
                            if 'free  energy   TOTEN' in line:
                                value = float(line.split('=')[-1].split()[0])
                                toten_values.append(value)
                    return toten_values[-1] if toten_values else None
                except Exception as e:
                    print(f"  ERROR: Failed to extract TOTEN in {folder_name}: {e}")
                    return None

            # Extract data
            length_c, z_lengths, plane_averages = read_z_locpot()
            e_fermi = read_fermi_energy()
            average_potential = calculate_average_potential(z_lengths, plane_averages)
            front_vacuum, _ = calculate_vacuum_averages(plane_averages)
            toten = extract_toten_from_outcar()

            # Skip if critical data is missing
            if None in [e_fermi, front_vacuum, average_potential, toten]:
                print(f"  Warning: Incomplete data in {folder_name}. Skipping.")
                os.chdir('..')
                continue

            # Append data
            data_list.append([e_fermi, front_vacuum, average_potential, toten, charge])

            os.chdir('..')

        except Exception as e:
            print(f"  ERROR: Failed to process {folder_name}: {e}")
            os.chdir('..')
            continue

    if not data_list:
        print("ERROR: No valid data extracted from any directories")
        print("Task 12 aborted.")
        return

    # Create and save DataFrame
    try:
        df = pd.DataFrame(data_list, columns=columns)
        output_file = os.path.join(task_output_dir, 'data_for_different_charges.xlsx')
        df.to_excel(output_file, index=False)
        print(f"\n✓ Task 12 completed successfully!")
        print(f"  Data saved to: {output_file}")
        print(f"  Processed {len(data_list)} directories")
    except Exception as e:
        print(f"ERROR: Failed to save Excel file: {e}")

    print("=" * 80)

def task13(args):
    """
    Task 13: ORR/OER Microkinetic Simulation and Visualization

    Performs microkinetic modeling for ORR/OER reactions, calculating coverage and current
    density as a function of potential (U vs. RHE). Generates data files and plots for
    O2 double-layer coverage, intermediate coverages, and polarization curves with key
    electrochemical parameters (half-wave potential, onset potential, etc.).

    Outputs:
    - U_coverage_current_data.txt: Coverage and current density data.
    - Key_potentials.txt: Key electrochemical parameters (U_1/2, U_10, U_onset).
    - O2_dl_Coverage.png: Plot of O2 double-layer coverage vs. potential.
    - Intermediate_Coverages.png: Plot of intermediate coverages vs. potential.
    - ORR_OER_Polarization_Curve.png: Polarization curve with annotated key potentials.
    """
    if not args.run_orr_oer_microkinetics:
        return

    print("=" * 80)
    print("TASK 13: ORR/OER MICROKINETIC SIMULATION".center(80))
    print("=" * 80)

    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_13_microkinetic_simulation"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    # --------------------------- Matplotlib Settings --------------------------- #
    plt.rcParams.update({
        'font.size': 14,
        'axes.labelweight': 'bold',
        'lines.linewidth': 2,
        'figure.autolayout': True
    })
    plt.switch_backend('agg')

    # ----------------------------- Constants from config ----------------------------- #
    config = MICROKINETICS_CONFIG
    acc = config["acc"]
    x_H2O = config["x_H2O"]
    x_O2 = config["x_O2"]
    T = config["T"]
    beta_i = np.array(config["beta_i"])
    e = config["e"]
    e_rho = config["e_rho"]
    kBT = config["kBT"]

    # Activation and reaction free energies from config
    delta_Ga0i = np.array(config["delta_Ga0i"])
    delta_Gi = np.array(config["delta_Gi"])
    Ai = np.array(config["Ai"])

    # File names
    RESULTS_FILE = os.path.join(task_output_dir, 'U_coverage_current_data.txt')
    KEY_FILE = os.path.join(task_output_dir, 'Key_potentials.txt')
    FIG_O2DL = os.path.join(task_output_dir, 'O2_dl_Coverage.png')
    FIG_COV = os.path.join(task_output_dir, 'Intermediate_Coverages.png')
    FIG_POLAR = os.path.join(task_output_dir, 'ORR_OER_Polarization_Curve.png')

    # U_range from config
    U_range = np.arange(config["U_range_start"], config["U_range_stop"], config["U_range_step"])

    # =========================== Utility Functions =========================== #

    def _K(idx: int, U: float) -> float:
        """Equilibrium constant K_i."""
        if idx <= 1:
            return float(np.exp(-delta_Gi[idx] / (kBT * T)))
        return float(np.exp(-(delta_Gi[idx] + U) / (kBT * T)))

    def _k_forward(idx: int, U: float) -> float:
        """Forward rate constant k_i."""
        if idx <= 1:
            return float(Ai[idx] * np.exp(-delta_Ga0i[idx] / (kBT * T)))
        return float(
            Ai[idx] * np.exp(-(delta_Ga0i[idx] + beta_i[idx] * (delta_Gi[idx] + U)) / (kBT * T))
        )

    def _k_backward(idx: int, U: float) -> float:
        """Backward rate constant k_-i."""
        K = _K(idx, U)
        kf = _k_forward(idx, U)
        return float(kf / K)

    def mkcrokinetic_equations(U: float) -> List[float]:
        """Solve coverage and calculate current density for a given potential U."""
        try:
            k1, k2, k3, k4, k5, k6 = (_k_forward(i, U) for i in range(6))
            k_1, k_2, k_3, k_4, k_5, k_6 = (_k_backward(i, U) for i in range(6))

            theta_O2_dl, theta_O2_star, theta_OOH_star, theta_O_star, theta_OH_star, theta_star = \
                symbols('theta_O2_dl theta_O2_star theta_OOH_star theta_O_star theta_OH_star theta_star')
            θ = (theta_O2_dl, theta_O2_star, theta_OOH_star, theta_O_star, theta_OH_star, theta_star)

            eq1 = k1 * x_O2 - k_1 * θ[0] - k2 * θ[0] * θ[5] + k_2 * θ[1]
            eq3 = N(-k3 * θ[1] + k_3 * θ[2] + k2 * θ[0] * θ[5] - k_2 * θ[1], acc)
            eq4 = N(-k4 * θ[2] + k_4 * θ[3] * x_H2O + k3 * θ[1] - k_3 * θ[2], acc)
            eq5 = N(-k5 * θ[3] + k_5 * θ[4] + k4 * θ[2] - k_4 * θ[3] * x_H2O, acc)
            eq6 = N(-k6 * θ[4] + k_6 * x_H2O * θ[5] + k5 * θ[3] - k_5 * θ[4], acc)
            eq7 = N(θ[5] + θ[4] + θ[3] + θ[2] + θ[1] - 1, acc)

            sol_all = solve((eq1, eq3, eq4, eq5, eq6, eq7), θ)
            sol = list(sol_all[-1])

            r1 = N(k1 * x_O2 - k_1 * sol[0], acc)
            r2 = N(k2 * sol[0] * sol[5] - k_2 * sol[1], acc)
            r3 = N(k3 * sol[1] - k_3 * sol[2], acc)
            r4 = N(k4 * sol[2] - k_4 * sol[3] * x_H2O, acc)

            current_density = - e_rho * (r1 + r2 + r3 + r4)

            sol.insert(0, U)
            sol.append(current_density)
            return sol
        except Exception as e:
            print(f"Error at U={U:.4f}V: {e}")
            return [U] + [0.0] * 6 + [0.0]

    def write_dat():
        """Generate and save simulation data."""
        results: List[List[float]] = []
        total_points = len(U_range)

        print(f"Starting ORR/OER microkinetic simulations for {total_points} potentials...")
        os.makedirs(args.output_dir, exist_ok=True)  # Ensure output directory exists

        for i, U in enumerate(U_range):
            row = mkcrokinetic_equations(U)
            results.append(row)
            if (i % 10 == 0) or (i == total_points - 1):
                print(f"Progress: {i + 1}/{total_points} - U={U:.3f}V, current={row[-1]:.3f} mA/cm²")

        with open(RESULTS_FILE, 'w') as f:
            header = (
                f"{'#U_vs_RHE':>12} {'theta_O2_dl':>12} {'theta_O2_star':>12} "
                f"{'theta_OOH_star':>12} {'theta_O_star':>12} {'theta_OH_star':>12} "
                f"{'theta_star':>12} {'current':>12}"
            )
            print(header, file=f)
            for row in results:
                print(" ".join(f"{v:12.6f}" for v in row), file=f)

        print(f"Simulation data saved to {RESULTS_FILE}")

    def _load_plot_dat(current_max: float = None) -> Tuple[np.ndarray, ...]:
        """Load simulation data and truncate at current_max."""
        if not os.path.exists(RESULTS_FILE):
            raise FileNotFoundError(f"Simulation data file {RESULTS_FILE} not found. Run write_dat() first.")

        try:
            data = np.loadtxt(RESULTS_FILE, comments="#")
            if current_max is not None:
                for i, current in enumerate(data[:, 7]):
                    if abs(current) > current_max:
                        data = data[:i]
                        break
            return (
                data[:, 0], data[:, 1], data[:, 2],
                data[:, 3], data[:, 4], data[:, 5],
                data[:, 6], data[:, 7]
            )
        except Exception as e:
            print(f"Error loading data: {e}")
            raise

    def plot_O_dl_coverage(current_max: float = None):
        """Plot O2 double-layer coverage vs. potential."""
        current_max = config["current_max"] if current_max is None else current_max
        try:
            U, theta_O2_dl, *_ = _load_plot_dat(current_max)
            plt.figure(figsize=(7, 6))
            smooth_U = np.linspace(U.min(), U.max(), 300)
            smooth_theta = interp1d(U, theta_O2_dl, kind='cubic')(smooth_U)

            plt.plot(smooth_U, smooth_theta, 'b-', linewidth=3, label='O₂ (dl)')
            plt.xlabel('U vs. RHE (V)')
            plt.ylabel('Coverage')
            plt.legend(loc='best')
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()
            plt.savefig(FIG_O2DL, dpi=500)
            plt.close()
            print(f"O2 dl coverage plot saved to {FIG_O2DL}")
        except Exception as e:
            print(f"Error plotting O2 dl coverage: {e}")

    def plot_intermediate_coverage(current_max: float = None):
        """Plot intermediate coverages vs. potential."""
        current_max = config["current_max"] if current_max is None else current_max
        try:
            U, _, theta_O2s, theta_OOH, theta_O, theta_OH, theta_star, _ = _load_plot_dat(current_max)
            plt.figure(figsize=(8, 6))

            plt.plot(U, theta_star, 'c-', label='*', linewidth=2)
            plt.plot(U, theta_O2s, 'm-', label='O₂*', linewidth=2)
            plt.plot(U, theta_OOH, 'g-', label='OOH*', linewidth=2)
            plt.plot(U, theta_O, 'r-', label='O*', linewidth=2)
            plt.plot(U, theta_OH, 'b-', label='OH*', linewidth=2)

            plt.xlabel('U vs. RHE (V)')
            plt.ylabel('Coverage')
            plt.legend(loc='best')
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()
            plt.savefig(FIG_COV, dpi=500)
            plt.close()
            print(f"Intermediate coverages plot saved to {FIG_COV}")
        except Exception as e:
            print(f"Error plotting intermediate coverages: {e}")

    def _find_nearest_idx(arr: np.ndarray, val: float) -> int:
        """Find index of nearest value in array."""
        return np.abs(arr - val).argmin()

    def plot_polarization_curve(current_max: float = None):
        """Plot polarization curve with key potentials."""
        current_max = config["current_max"] if current_max is None else current_max
        try:
            U, *_, current = _load_plot_dat(current_max)

            half_wave_current = current[0] / 2
            current_10mA = 10.0
            onset_current = -0.01

            U_half = U[_find_nearest_idx(current, half_wave_current)]
            U_10 = U[_find_nearest_idx(current, current_10mA)]
            U_onset = U[_find_nearest_idx(current, onset_current)]

            print("\n=== Key Electrochemical Parameters ===")
            print(f"Half-wave potential (U_½) = {U_half:.3f} V")
            print(f"Potential at 10 mA/cm² (U_10) = {U_10:.3f} V")
            print(f"Onset potential (U_onset) = {U_onset:.3f} V")

            with open(KEY_FILE, 'w') as f:
                f.write(f"U_1/2 = {U_half:.3f} V\n")
                f.write(f"U_10 = {U_10:.3f} V\n")
                f.write(f"U_onset = {U_onset:.3f} V\n")

            plt.figure(figsize=(9, 7))
            plt.plot(U, current, 'b-', linewidth=3, label='ORR/OER Current')
            plt.axhline(y=0, color='k', linestyle='-', alpha=0.5)

            def add_annotation(x, y, label, color):
                plt.axhline(y=y, color=color, linestyle='--', alpha=0.7)
                plt.axvline(x=x, color=color, linestyle='--', alpha=0.7)
                plt.annotate(
                    f'{label} = {x:.2f} V',
                    xy=(x, y),
                    xytext=(x + 0.02, y + 2),
                    arrowprops=dict(arrowstyle='->', color=color),
                    fontsize=12,
                    color=color
                )

            add_annotation(U_half, half_wave_current, r'$U_{1/2}$', 'r')
            add_annotation(U_10, current_10mA, r'$U_{10}$', 'g')
            add_annotation(U_onset, onset_current, r'$U_{onset}$', 'b')

            if np.any(current < 0):
                plt.text(0.8, -5, 'ORR Region', fontsize=14, ha='center',
                         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="blue", lw=1))
            if np.any(current > 0):
                plt.text(1.1, 5, 'OER Region', fontsize=14, ha='center',
                         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red", lw=1))

            plt.xlabel('U vs. RHE (V)')
            plt.ylabel('Current Density (mA/cm²)')
            plt.title('ORR/OER Polarization Curve')
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(FIG_POLAR, dpi=500)
            plt.close()
            print(f"Polarization curve plot saved to {FIG_POLAR}")
        except Exception as e:
            print(f"Error plotting polarization curve: {e}")

    # Main execution
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        write_dat()
        print("\nGenerating plots...")
        plot_O_dl_coverage()
        plot_intermediate_coverage()
        plot_polarization_curve()
        print("\n✓ Task 13 completed successfully!")
        print(f"Results saved to: {args.output_dir}")
    except Exception as e:
        print(f"ERROR: Task 13 failed: {e}")

    print("=" * 80)





def task14(args):
    # ✅ Allow Task 14 to run if explicitly selected via --task 14
    if not args.run_oer_orr and 14 not in args.task:
        print("Task 14 skipped: neither --run-oer-orr nor --task 14 specified")
        return


    print("=" * 80)
    print("TASK 14: OER/ORR Free Energy Diagram (Constant Charge)".center(80))
    print("=" * 80)
    task_output_dir = os.path.join(
        args.output_dir,
        "TASK_14_constant_charge_diagram"
    )
    os.makedirs(task_output_dir, exist_ok=True)

    # 创建子目录
    oer_orr_dir = create_output_directory(os.path.join(task_output_dir, "oer_orr_analysis"))
    print(f"Created directory: {oer_orr_dir}")

    # 加载配置
    config = OER_ORR_CONFIG.copy()
    config["pH"] = args.pH if args.pH is not None else config["pH"]
    config["reaction_type"] = args.reaction_type if args.reaction_type else config["reaction_type"]
    pH = config["pH"]
    U = args.U if args.U is not None else config.get("U", 0.0)  # 默认电位为 0 V
    print(f"Config: pH={pH}, U={U} V, reaction_type={config['reaction_type']}")

    # 读取 VASP 数据和矫正值
    energies = config["vasp_energies"]
    corrections = config["corrections"]
    kT = 0.0257  # eV at 298 K
    e = 1.0      # eV/V

    E_slab = energies["E_slab"]
    E_OOH = energies["E_OOH"]
    E_O = energies["E_O"]
    E_OH = energies["E_OH"]
    E_H2O = energies["E_H2O"]
    E_O2 = energies["E_O2"]
    E_H2 = energies["E_H2"]
    delta_G_water = corrections["delta_G_water"]
    ZPE_TDS = corrections["ZPE_TDS"]

    # 计算吸附自由能 (基础值，无 pH/U 修正)
    delta_GOOH = (E_OOH - E_slab - (2 * E_H2O - 1.5 * E_H2)) + ZPE_TDS["*OOH"]
    delta_GO = (E_O - E_slab - (E_H2O - E_H2)) + ZPE_TDS["*O"]
    delta_GOH = (E_OH - E_slab - (E_H2O - 0.5 * E_H2)) + ZPE_TDS["*OH"]
    print(f"Base intermediates: ΔG*OOH={delta_GOOH:.4f} eV, ΔG*O={delta_GO:.4f} eV, ΔG*OH={delta_GOH:.4f} eV")

    # 计算台阶能，加入 pH 和 U 修正
    if config["reaction_type"] == "ORR":
        dG_steps = [
            delta_GOOH  - delta_G_water - (-kT * np.log(10) * pH) + (e * U),  # * + O2 + H+ + e- → *OOH
            (delta_GO - delta_GOOH) - (-kT * np.log(10) * pH) + (e * U),  # *OOH + H+ + e- → *O + H2O
            (delta_GOH - delta_GO) - (-kT * np.log(10) * pH) + (e * U),   # *O + H+ + e- → *OH
            (-delta_GOH ) - (-kT * np.log(10) * pH) + (e * U)  # *OH + H+ + e- → * + H2O (adjusted)
        ]
        initial_height = delta_G_water - 4 * (kT * np.log(10) * pH + e * U)  # 从水形成的总自由能减去 4e- 贡献
    else:  # OER
        dG_steps = [
            delta_GOH + (kT * np.log(10) * pH) - (e * U),               # * + H2O → *OH + H+ + e- (逆向)
            (delta_GO - delta_GOH) + (kT * np.log(10) * pH) - (e * U),  # *OH → *O + H+ + e- (逆向)
            (delta_GOOH - delta_GO) + (kT * np.log(10) * pH) - (e * U), # *O + H2O → *OOH + H+ + e- (逆向)
            (delta_G_water - delta_GOOH) + (kT * np.log(10) * pH) - (e * U)  # *OOH → * + O2 + H+ + e- (逆向)
        ]
        initial_height = 0  # OER 从基态开始

    print(f"Step energies with pH/U correction: ΔG1={dG_steps[0]:.4f} eV, ΔG2={dG_steps[1]:.4f} eV, "
          f"ΔG3={dG_steps[2]:.4f} eV, ΔG4={dG_steps[3]:.4f} eV")
    print(f"Initial height for {config['reaction_type']}: {initial_height:.4f} eV")

    # 累积自由能
    cumulative_energy = [initial_height]
    current_energy = initial_height
    for dG in dG_steps:
        current_energy += dG
        cumulative_energy.append(current_energy)

    # 计算反应势垒
    barriers = []
    for i in range(4):
        barrier = cumulative_energy[i + 1] - cumulative_energy[i]
        barriers.append(barrier)

    rds_idx = np.argmax(barriers)
    step_names = ["Step 1", "Step 2", "Step 3", "Step 4"]
    step_desc = {
        "ORR": ["* + O2 + H+ + e- → *OOH", "*OOH + H+ + e- → *O + H2O", "*O + H+ + e- → *OH", "*OH + H+ + e- → * + H2O"],
        "OER": ["* + H2O → *OH + H+ + e-", "*OH → *O + H+ + e-", "*O + H2O → *OOH + H+ + e-", "*OOH → * + O2 + H+ + e-"]
    }[config["reaction_type"]]
    rds_step_name = step_names[rds_idx]
    rds_step_desc = step_desc[rds_idx]
    rds_barrier = barriers[rds_idx]
    rds_end_energy = cumulative_energy[rds_idx + 1]

    # 绘制路径图
    plt.figure(figsize=(10, 8))
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, 1))  # 单条件

    labels = [
        ["* + O₂", "*OOH", "*O", "*OH", "* + H₂O"],
        ["* + H2O", "*OH", "*O", "*OOH", "* + O2"]
    ][0 if config["reaction_type"] == "ORR" else 1]

    x_points = []
    y_points = []
    for i in range(5):
        x_plat_start = i - 0.4
        x_plat_end = i + 0.4
        y_plat = cumulative_energy[i]
        x_points.extend([x_plat_start, x_plat_end])
        y_points.extend([y_plat, y_plat])

        if i < 4:
            x_slope_start = x_plat_end
            x_slope_end = i + 1 - 0.4
            x_mid = (x_slope_start + x_slope_end) / 2
            x_points.extend([x_slope_start, x_mid])
            y_points.extend([y_plat, cumulative_energy[i + 1]])
            x_points.extend([x_mid, x_slope_end])
            y_points.extend([cumulative_energy[i + 1], cumulative_energy[i + 1]])

    plt.plot(
        x_points, y_points,
        color=colors[0],
        linewidth=2.5,
        label=f"RDS: {rds_step_name}, Barrier: {rds_barrier:.2f} eV"
    )

    for i, energy in enumerate(cumulative_energy):
        plt.text(i, energy + 0.05, f"{energy:.2f}", ha='center', fontsize=9, color=colors[0])

    plt.xticks(range(5), labels, fontsize=12)
    plt.xlim(-0.5, 4.5)
    plt.ylabel('Free Energy (eV)', fontsize=14)
    plt.title(f'{config["reaction_type"]} Reaction Pathway (pH={pH}, U={U} V)', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10, loc='best', framealpha=0.7)
    plt.tight_layout()

    fig_path = os.path.join(oer_orr_dir, f"{config['reaction_type']}_pathway_pH{pH}_U{U}.pdf")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.savefig(fig_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Pathway diagram saved: {fig_path}")
    plt.close()

    # 保存数据
    intermediates_df = pd.DataFrame({"delta_GOOH": [delta_GOOH], "delta_GO": [delta_GO], "delta_GOH": [delta_GOH]}, index=[U])
    steps_df = pd.DataFrame({"d_G1": [dG_steps[0]], "d_G2": [dG_steps[1]], "d_G3": [dG_steps[2]], "d_G4": [dG_steps[3]]}, index=[U])
    intermediates_path = os.path.join(oer_orr_dir, f'intermediates_pH{pH}_{config["reaction_type"]}_data.csv')
    steps_path = os.path.join(oer_orr_dir, f'steps_pH{pH}_{config["reaction_type"]}_data.csv')
    intermediates_df.to_csv(intermediates_path)
    steps_df.to_csv(steps_path)
    print(f"Intermediates data saved to {intermediates_path}")
    print(f"Steps data saved to {steps_path}")

    # 保存 RDS 信息
    rds_path = os.path.join(oer_orr_dir, f"{config['reaction_type']}_rds_info_pH{pH}_U{U}.csv")
    with open(rds_path, 'w', encoding='utf-8') as f:
        f.write("pH,U,RDS_Step,RDS_Description,RDS_Barrier(eV),RDS_End_Energy(eV)\n")
        f.write(f"{pH},{U},{rds_step_name},\"{rds_step_desc}\",{rds_barrier},{rds_end_energy}\n")
    print(f"RDS details saved: {rds_path}")

    print(f"\n✓ {config['reaction_type']} pathway analysis completed. Results saved to: {oer_orr_dir}")


def task15(args):
    """
    Task 15: OER/ORR Activity Volcano Plot
    Optimized boundary line handling with full data export
    Generates volcano plots for both OER and ORR reactions,
    saving results to separate subdirectories under '2D_Volcano'.
    Incorporates inverse relationship between OER and ORR free energy steps.
    """
    if not args.run_oer_orr_volcano:
        return

    print("=" * 80)
    print("TASK 15: OER/ORR ACTIVITY VOLCANO PLOT".center(80))
    print("=" * 80)

    # 统一处理 OER 和 ORR
    for reaction_type in ["OER", "ORR"]:
        print(f"\n--- Generating {reaction_type.upper()} Volcano Plot ---")

        # 创建输出目录
        task15_output_dir = os.path.join(
            args.output_dir,
            "TASK_15_oer_orr_volcano"
        )
        base_output_dir = create_output_directory(
            os.path.join(task15_output_dir, "2D_Volcano")
        )

        current_reaction_output_dir = create_output_directory(os.path.join(base_output_dir, reaction_type.upper()))

        # 提取金属吸附能数据
        metal_data = {}
        for metal, energies in METAL_ADSORPTION_ENERGIES.items():
            delta_G_OH, delta_G_O, delta_G_OOH = energies
            metal_data[metal] = {
                'delta_G_OH': delta_G_OH,
                'delta_G_O': delta_G_O,
                'delta_G_OOH': delta_G_OOH,
                'delta_G_O_minus_OH': delta_G_O - delta_G_OH
            }

        # ΔG_OOH vs ΔG_OH 线性拟合 (此拟合关系不随反应类型改变)
        delta_G_OH_vals = [d['delta_G_OH'] for d in metal_data.values()]
        delta_G_OOH_vals = [d['delta_G_OOH'] for d in metal_data.values()]
        # 检查是否所有 delta_G_OH_vals 都相同，避免 np.polyfit 警告
        if len(set(delta_G_OH_vals)) == 1:
            slope_OOH = 0
            intercept_OOH = delta_G_OOH_vals[0] if delta_G_OOH_vals else 0
            print(
                "Warning: All delta_G_OH_vals are identical, cannot perform linear fit for delta_G_OOH. Using constant value.")
        else:
            slope_OOH, intercept_OOH = np.polyfit(delta_G_OH_vals, delta_G_OOH_vals, 1)  # m, c

        # 设置坐标范围
        if VOLCANO_PLOT_RANGE:
            x_min, x_max, y_min, y_max = VOLCANO_PLOT_RANGE
            print(f"Using configured plot range: x=[{x_min},{x_max}], y=[{y_min},{y_max}]")
        else:
            margin = 0.5
            x_min = min(delta_G_OH_vals) - margin
            x_max = max(delta_G_OH_vals) + margin
            y_min = min([d['delta_G_O_minus_OH'] for d in metal_data.values()]) - margin
            y_max = max([d['delta_G_O_minus_OH'] for d in metal_data.values()]) + margin
            print(f"Using auto-calculated plot range: x=[{x_min},{x_max}], y=[{y_min},{y_max}]")

        resolution = 300
        x = np.linspace(x_min, x_max, resolution)
        y = np.linspace(y_min, y_max, resolution)
        X, Y = np.meshgrid(x, y)

        # 常数定义
        # total_g 代表 ΔG_total_4e, 即 4.92 eV
        total_g = 4.92
        equilibrium_potential = 1.23  # 理论平衡电位

        # 计算自由能步骤 (基于网格 X=dG_OH_mesh, Y=dG_O_minus_OH_mesh)
        dG_OH_mesh = X
        dG_O_minus_OH_mesh = Y
        dG_O_mesh = dG_OH_mesh + dG_O_minus_OH_mesh  # ΔG_O = ΔG_OH + (ΔG_O - ΔG_OH)
        dG_OOH_mesh = slope_OOH * dG_OH_mesh + intercept_OOH  # ΔG_OOH 由线性拟合得到

        # 根据反应类型定义自由能步骤和计算过电位
        if reaction_type == "OER":
            # OER 自由能步骤
            d_G1 = dG_OH_mesh  # H2O -> OH*
            d_G2 = dG_O_minus_OH_mesh  # OH* -> O*
            d_G3 = dG_OOH_mesh - dG_O_mesh  # O* -> OOH*
            d_G4 = total_g - dG_OOH_mesh  # OOH* -> O2 + H2O

            max_dG = np.maximum.reduce([d_G1, d_G2, d_G3, d_G4])
            overpotentials = np.maximum(max_dG - equilibrium_potential, 0)
            title_suffix = "OER"
            cbar_label_suffix = r'Overpotential $\eta_{OER}$ (V)'
            plot_filename_prefix = "OER_activity_volcano"

        elif reaction_type == "ORR":
            # ORR 自由能步骤 (根据你的建议，与OER步骤呈逆关系)
            d_G1 = dG_OOH_mesh - total_g  # O2 + H2O -> OOH* (对应 -dG4_OER)
            d_G2 = dG_O_mesh - dG_OOH_mesh  # OOH* -> O* (对应 -dG3_OER)
            d_G3 = dG_OH_mesh - dG_O_mesh  # O* -> OH* (对应 -dG2_OER)
            d_G4 = -dG_OH_mesh  # OH* -> H2O (对应 -dG1_OER, 假设G_H2O=0)

            max_dG = np.maximum.reduce([d_G1, d_G2, d_G3, d_G4])
            # ORR的过电位，同样采用 max(ΔG_i) - E_eq 形式，使得火山图形状一致
            overpotentials = np.maximum(max_dG + equilibrium_potential,
                                        0)  # 修正：ORR通常是max(dG_i - U_0)或者max(abs(dG_i - U_0))
            # 如果是dG_i是负值，max_dG + equilibrium_potential会是负的。
            # 为了形成火山图，我们通常希望中心是低的。
            # 假设d_G1, d_G2, d_G3, d_G4 已经做了调整，是正向的能量垒
            # 否则，这里需要重新考虑 ORR overpotential 的计算。
            # 当前你提供的 d_G 定义，max_dG + equilibrium_potential 应该能得到正值。

            title_suffix = "ORR"
            cbar_label_suffix = r'Overpotential $\eta_{ORR}$ (V)'
            plot_filename_prefix = "ORR_activity_volcano"

        # 绘图
        fig, ax = plt.subplots(figsize=(8, 6))

        # 设置颜色范围
        if OVERPTENTIAL_COLOR_RANGE:
            vmin, vmax = OVERPTENTIAL_COLOR_RANGE
            print(f"Using configured color range: vmin={vmin}, vmax={vmax}")
        else:
            # 自动调整 vmax 到过电位数据的 95th percentile
            if np.any(overpotentials > 0):  # 确保有过电位值
                vmax = np.percentile(overpotentials[overpotentials > 0], 95) * 1.2  # 稍微扩大范围
                vmin = np.min(overpotentials[overpotentials > 0])  # 或者设置一个最小值
            else:
                vmin, vmax = 0, 1.5  # 默认值
            print(f"Adjusted color range: vmin={vmin}, vmax={vmax}")

        n_levels = 500
        levels = np.linspace(vmin, vmax, n_levels)

        # 绘制火山图
        contour = ax.contourf(X, Y, overpotentials,
                              levels=levels,
                              cmap='hot',  # 使用逆向颜色映射，使得低过电位（蓝色）更突出
                              alpha=0.8,
                              vmin=vmin,
                              vmax=vmax)
        cbar = plt.colorbar(contour, ax=ax)
        cbar.set_label(cbar_label_suffix)

        # 添加金属数据点
        for metal, data in metal_data.items():
            x_val = data['delta_G_OH']
            y_val = data['delta_G_O_minus_OH']
            ax.scatter(x_val, y_val, c='cyan', edgecolors='black',
                       s=120, zorder=3, alpha=0.9)
            ax.text(x_val, y_val + 0.05, metal,
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

        # ===== 关键边界线处理 =====
        m = slope_OOH
        c = intercept_OOH

        # 存储边界线数据用于导出
        boundary_line_points = []

        # 根据反应类型定义边界线
        if reaction_type == "OER":
            lines = [
                {'pair': (1, 2), 'slope': 1.0, 'intercept': 0.0, 'gi_index': 1, 'conditions': [3, 4]},  # dG1 = dG2
                {'pair': (2, 3), 'slope': (m - 1) / 2, 'intercept': c / 2, 'gi_index': 2, 'conditions': [1, 4]},
                # dG2 = dG3
                {'pair': (3, 4), 'slope': 2 * m - 1, 'intercept': 2 * c - total_g, 'gi_index': 4, 'conditions': [1, 2]},
                # dG3 = dG4
                # 非相邻边界线
                {'pair': (1, 3), 'slope': m - 2, 'intercept': c, 'gi_index': 1, 'conditions': [2, 4]},  # dG1 = dG3
                {'pair': (2, 4), 'slope': -m, 'intercept': total_g - c, 'gi_index': 2, 'conditions': [1, 3]},
                # dG2 = dG4
            ]
        elif reaction_type == "ORR":
            # 根据你提供的ORR步骤定义 d_G1, d_G2, d_G3, d_G4 重新推导边界线
            # ORR 自由能步骤 (这里引用的是上面 elif reaction_type == "ORR": 中定义的 d_G1, ..., d_G4)
            # G1_orr = dG_OOH_mesh - total_g
            # G2_orr = dG_O_mesh - dG_OOH_mesh
            # G3_orr = dG_OH_mesh - dG_O_mesh
            # G4_orr = -dG_OH_mesh

            lines = [
                # 1. dG1 = dG2
                {'pair': (1, 2), 'slope': (2 * m - 1), 'intercept': (2 * c - total_g), 'gi_index': 1,
                 'conditions': [3, 4]},
                # 2. dG2 = dG3
                {'pair': (2, 3), 'slope': (m - 1) / 2, 'intercept': c / 2, 'gi_index': 2, 'conditions': [1, 4]},
                # 3. dG3 = dG4
                {'pair': (3, 4), 'slope': 1.0, 'intercept': 0.0, 'gi_index': 3, 'conditions': [1, 2]},
                # 4. G4_orr = G1_orr (垂直线，特殊处理)
                #    -dG_OH_mesh = dG_OOH_mesh - total_g => -X = mX + c - total_g => X = (total_g - c) / (1 + m)
                #    这里我们不添加到 lines 列表，而是在循环外单独绘制
                # 5. dG1 = dG3
                {'pair': (1, 3), 'slope': -m, 'intercept': total_g - c, 'gi_index': 1, 'conditions': [2, 4]},
                # 6. dG2 = dG4
                {'pair': (2, 4), 'slope': (m - 2), 'intercept': c, 'gi_index': 2, 'conditions': [1, 3]},
            ]

        # 绘制边界线
        for line in lines:
            slope = line['slope']
            intercept = line['intercept']
            gi_index = line['gi_index']
            cond_k, cond_l = line['conditions']

            x_line = np.linspace(x_min, x_max, 3000)
            y_line = slope * x_line + intercept

            mask_clip = (y_line >= y_min) & (y_line <= y_max)
            x_line_c = x_line[mask_clip]
            y_line_c = y_line[mask_clip]

            if len(x_line_c) == 0:
                continue

            # 根据当前 reaction_type 重新计算沿线的自由能步骤
            if reaction_type == "OER":
                G1_line = x_line_c
                G2_line = y_line_c
                G3_line = (m * x_line_c + c) - (x_line_c + y_line_c)
                G4_line = total_g - (m * x_line_c + c)
            elif reaction_type == "ORR":
                # 确保这里使用的 G_line 定义与 ORR overpotentials 的定义一致
                G1_line = (m * x_line_c + c) - total_g
                G2_line = (x_line_c + y_line_c) - (m * x_line_c + c)
                G3_line = x_line_c - (x_line_c + y_line_c)
                G4_line = -x_line_c

            all_G_lines_at_boundary = [G1_line, G2_line, G3_line, G4_line]

            Gi_line_ref = all_G_lines_at_boundary[gi_index - 1]

            # 检查在当前线段上，参考步骤是否是最大的（与cond_k, cond_l定义的步骤比较）
            mask_boundary = (Gi_line_ref >= all_G_lines_at_boundary[cond_k - 1]) & \
                            (Gi_line_ref >= all_G_lines_at_boundary[cond_l - 1])

            x_seg = x_line_c[mask_boundary]
            y_seg = y_line_c[mask_boundary]

            if len(x_seg) == 0:
                continue

            # 找到连续的线段并绘制
            diff_mask = np.diff(mask_boundary.astype(int))
            starts = np.where(diff_mask == 1)[0] + 1
            ends = np.where(diff_mask == -1)[0] + 1

            if mask_boundary[0]:
                starts = np.insert(starts, 0, 0)
            if mask_boundary[-1]:
                ends = np.append(ends, len(x_line_c))

            if len(starts) > len(ends):
                ends = np.append(ends, len(x_line_c))
            elif len(ends) > len(starts):
                starts = np.insert(starts, 0, 0)

            for s, e in zip(starts, ends):
                if s < e:
                    ax.plot(x_line_c[s:e], y_line_c[s:e], 'k--', lw=1.5, alpha=0.8)
                    if len(x_line_c[s:e]) > 1:
                        boundary_line_points.append({
                            'pair': f"{line['pair'][0]}-{line['pair'][1]}",
                            'slope': slope,
                            'intercept': intercept,
                            'x_start': x_line_c[s],
                            'y_start': y_line_c[s],
                            'x_end': x_line_c[e - 1],
                            'y_end': y_line_c[e - 1]
                        })

        # 特殊处理 ORR 的垂直边界线 G4_orr = G1_orr
        if reaction_type == "ORR":
            if 1 + m != 0:  # 避免除以零
                x_vertical = (total_g - c) / (1 + m)
                if x_min <= x_vertical <= x_max:
                    # 确定这条线段的有效Y范围，它也必须处于限速区域
                    # G1_orr = G4_orr 意味着 Gi_line_ref 实际上就是 G1_orr 或 G4_orr
                    # 我们需要检查 G1_orr >= G2_orr 且 G1_orr >= G3_orr

                    # 重新计算沿这条垂直线的 G1, G2, G3, G4
                    # 这里的 x_line_c 只有一个值 x_vertical
                    # Y 对应于 y_line，所以我们需要一个 y_range
                    y_vert_line = np.linspace(y_min, y_max, 300)

                    G1_vert_line = (m * x_vertical + c) - total_g
                    G2_vert_line = (x_vertical + y_vert_line) - (m * x_vertical + c)
                    G3_vert_line = x_vertical - (x_vertical + y_vert_line)
                    G4_vert_line = -x_vertical

                    all_G_vert_lines = [G1_vert_line, G2_vert_line, G3_vert_line, G4_vert_line]

                    # 这里的 ref_step_idx 可以是 0 (G1) 或 3 (G4)，因为我们检查的是 G1=G4
                    # 条件是 G1 >= G2 且 G1 >= G3
                    mask_boundary_vert = (G1_vert_line >= all_G_vert_lines[1]) & \
                                         (G1_vert_line >= all_G_vert_lines[2])

                    y_seg_vert = y_vert_line[mask_boundary_vert]

                    if len(y_seg_vert) > 0:
                        # 找到连续的线段并绘制
                        diff_mask_vert = np.diff(mask_boundary_vert.astype(int))
                        starts_vert = np.where(diff_mask_vert == 1)[0] + 1
                        ends_vert = np.where(diff_mask_vert == -1)[0] + 1

                        if mask_boundary_vert[0]:
                            starts_vert = np.insert(starts_vert, 0, 0)
                        if mask_boundary_vert[-1]:
                            ends_vert = np.append(ends_vert, len(y_vert_line))

                        if len(starts_vert) > len(ends_vert):
                            ends_vert = np.append(ends_vert, len(y_vert_line))
                        elif len(ends_vert) > len(starts_vert):
                            starts_vert = np.insert(starts_vert, 0, 0)

                        for s_v, e_v in zip(starts_vert, ends_vert):
                            if s_v < e_v:
                                ax.plot([x_vertical] * len(y_vert_line[s_v:e_v]), y_vert_line[s_v:e_v], 'k--', lw=1.5,
                                        alpha=0.8)
                                if len(y_vert_line[s_v:e_v]) > 1:
                                    boundary_line_points.append({
                                        'pair': '4-1 (vertical)',
                                        'slope': np.inf,  # 无穷斜率表示垂直线
                                        'intercept': x_vertical,  # 截距在此表示 x 值
                                        'x_start': x_vertical,
                                        'y_start': y_vert_line[s_v],
                                        'x_end': x_vertical,
                                        'y_end': y_vert_line[e_v - 1]
                                    })
            else:
                print("Warning: For ORR, 1 + m is zero, vertical boundary line G4=G1 cannot be calculated.")

        # 设置图形标题和标签
        ax.set_title(f'{title_suffix} Activity Volcano Map', fontsize=14)
        ax.set_xlabel(r'$\Delta G_{OH^*}$ (eV)')
        ax.set_ylabel(r'$\Delta G_{O^*} - \Delta G_{OH^*}$ (eV)')

        # 设置坐标轴范围
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        plt.tight_layout()

        # 保存图表
        plot_filename = os.path.join(current_reaction_output_dir, f"{plot_filename_prefix}.png")
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved {reaction_type.upper()} volcano plot to: {plot_filename}")

        # ===== 数据导出部分 =====
        # 1. 导出火山网格数据
        X_flat = X.flatten()
        Y_flat = Y.flatten()
        overpotentials_flat = overpotentials.flatten()

        grid_df = pd.DataFrame({
            'x_delta_G_OH': X_flat,
            'y_delta_G_O_minus_OH': Y_flat,
            'overpotential': overpotentials_flat
        })
        grid_filename = os.path.join(current_reaction_output_dir, f"{reaction_type.lower()}_volcano_grid_data.csv")
        grid_df.to_csv(grid_filename, index=False)
        print(f"Saved grid data to: {grid_filename}")

        # 2. 导出金属点数据
        metal_list = []
        for metal, data in metal_data.items():
            # 计算金属点的过电位，同样使用修正后的 ΔG 步骤
            dG_OH_val = data['delta_G_OH']
            dG_O_minus_OH_val = data['delta_G_O_minus_OH']
            dG_O_val = dG_OH_val + dG_O_minus_OH_val
            dG_OOH_val = slope_OOH * dG_OH_val + intercept_OOH

            if reaction_type == "OER":
                d_G1_metal = dG_OH_val
                d_G2_metal = dG_O_minus_OH_val
                d_G3_metal = dG_OOH_val - dG_O_val
                d_G4_metal = total_g - dG_OOH_val
            elif reaction_type == "ORR":
                d_G1_metal = dG_OOH_val - total_g
                d_G2_metal = dG_O_val - dG_OOH_val
                d_G3_metal = dG_OH_val - dG_O_val
                d_G4_metal = -dG_OH_val

            max_step_metal = max(d_G1_metal, d_G2_metal, d_G3_metal, d_G4_metal)
            # 确保金属点过电位计算与网格过电位计算一致
            eta = np.maximum(max_step_metal + equilibrium_potential, 0)  # 修正：这里也要使用 + equilibrium_potential

            metal_list.append({
                'metal': metal,
                'x_delta_G_OH': dG_OH_val,
                'y_delta_G_O_minus_OH': dG_O_minus_OH_val,
                'eta_overpotential': eta
            })

        metal_df = pd.DataFrame(metal_list)
        metal_filename = os.path.join(current_reaction_output_dir, f"{reaction_type.lower()}_volcano_metal_points.csv")
        metal_df.to_csv(metal_filename, index=False)
        print(f"Saved metal point data to: {metal_filename}")

        # 3. 导出边界线数据
        if boundary_line_points:
            boundaries_df = pd.DataFrame(boundary_line_points)
            boundary_filename = os.path.join(current_reaction_output_dir,
                                             f"{reaction_type.lower()}_volcano_boundary_lines.csv")
            boundaries_df.to_csv(boundary_filename, index=False)
            print(f"Saved boundary line data to: {boundary_filename}")
        else:
            print("Warning: No boundary line data to export for " + reaction_type.upper())

        # 4. 导出线性拟合参数
        fit_data = pd.DataFrame({
            'parameter': ['slope_OOH', 'intercept_OOH'],
            'value': [slope_OOH, intercept_OOH]
        })
        fit_filename = os.path.join(current_reaction_output_dir, f"{reaction_type.lower()}_volcano_linear_fit.csv")
        fit_data.to_csv(fit_filename, index=False)
        print(f"Saved linear fit parameters to: {fit_filename}")

        print(f"--- {reaction_type.upper()} Volcano Plot Generation Complete ---")

    print("\n" + "-" * 80)
    print("TASK 15 COMPLETED SUCCESSFULLY")
    print("-" * 80)

# Update TASK_REGISTRY
TASK_REGISTRY = {
    1: task1,
    2: task2,
    3: task3,
    4: task4,
    5: task5,
    6: task6,
    7: task7,
    8: task8,
    9: task9,
    10: task10,
    11: task11,
    12: task12,
    13: task13,
    14: task14,
    15:task15
      # New task for Pourbaix diagram
}
