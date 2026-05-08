import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# 包内导入
from cp_analyzer.utils import safe_filename, create_output_directory, E_free_vs_SHE, quadratic_function_vs_she, U_rhe_to_U_she
from scipy.optimize import curve_fit
from cp_analyzer.config_default import LITERATURE_PARAMS
import os
from cp_analyzer.data_processing import Q_tran_vs_pH, read_excel_columns_by_index, calculate_ads_energy

def plot_save_Q_vs_pH(output_dir, file_paths=None, Urhe=0, degree=3):
    """Plot and save Q vs pH relation (two plots)"""
    if not file_paths:
        print("Warning: No file paths provided, skipping dQ vs pH analysis")
        return

    Q, pH = Q_tran_vs_pH(file_paths=file_paths, Urhe=Urhe)

    if not Q or not pH or len(Q) == 0 or len(pH) == 0:
        print("Warning: No valid dQ/pH data, skipping dQ vs pH analysis,you need input dQ(the charge change) \n of "
              " intermediates in *.xlsx at the last column!!")
        return

    # Determine system names (based on whether O2 is included)
    include_O2 = any('slab.xlsx' in p.lower() for p in file_paths)
    if include_O2 and len(file_paths) >= 5:
        system_names = {0: 'Slab', 1: 'O2', 2: 'OOH', 3: 'O', 4: 'OH'}
    elif len(file_paths) >= 4:
        system_names = {0: 'Slab', 1: 'OOH', 2: 'O', 3: 'OH'}
    else:
        system_names = {i: f"System {i + 1}" for i in range(len(file_paths))}

    fit_params = []
    fit_accuracies = []

    # Create output directory
    plot_dir = create_output_directory(os.path.join(output_dir, "Q_vs_pH"))

    # First plot: Q vs pH fitted curves
    plt.figure(figsize=(10, 7))
    for i, (pH_values, Q_values) in enumerate(zip(pH, Q)):
        valid_mask = np.isfinite(pH_values) & np.isfinite(Q_values)
        pH_values = pH_values[valid_mask]
        Q_values = Q_values[valid_mask]

        if len(pH_values) < degree + 1:
            print(f"Warning: Not enough points for {system_names.get(i, 'System')} {i + 1}")
            continue

        try:
            poly_params = np.polyfit(pH_values, Q_values, degree)
            Q_fit = np.polyval(poly_params, pH_values)
            ss_res = np.sum((Q_values - Q_fit) ** 2)
            ss_tot = np.sum((Q_values - np.mean(Q_values)) ** 2) if np.std(Q_values) > 0 else 1.0
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        except Exception as e:
            print(f"Polynomial fit error: {e}")
            continue

        fit_params.append(poly_params)
        fit_accuracies.append(r_squared)

        # Calculate fitted curve
        pH_fit = np.linspace(min(pH_values), max(pH_values), 100)
        Q_fit = np.polyval(poly_params, pH_fit)

        # Plot fitted curve and original data points for current system
        plt.plot(pH_fit, Q_fit, label=f'{system_names.get(i, f"System {i + 1}")} Fit', linewidth=2)
        plt.scatter(pH_values, Q_values, s=60, alpha=0.7)

        # Save fitted data
        fit_data = np.column_stack((pH_fit, Q_fit))
        system_label = system_names.get(i, f"System{i + 1}").replace('*', '').replace(' ', '_')
        fit_file = os.path.join(plot_dir, f'fit_data_{system_label}.txt')
        np.savetxt(fit_file, fit_data, header='pH, Q', fmt='%.6f')
        print(f"Saved fit data to: {fit_file}")

    plt.xlabel('pH', fontsize=12)
    plt.ylabel('Q', fontsize=12)
    plt.title('Q vs pH Polynomial Fit', fontsize=14)
    plt.legend(fontsize=10, loc='best')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    fit_plot = os.path.join(plot_dir, 'Q_vs_pH_fit.png')
    plt.savefig(fit_plot, dpi=300)
    plt.close()
    print(f"Saved Q vs pH fit plot: {fit_plot}")

    # Second plot: Q values in pH 0-14 range
    if not fit_params:
        print("Warning: No valid polynomial fits found, skipping pH range plot")
        return

    plt.figure(figsize=(10, 7))
    pH_range_values = np.linspace(0, 14, 15)  # Create pH range 0-14
    Q_pH_range = np.column_stack([np.polyval(params, pH_range_values) for params in fit_params])

    for col_idx in range(Q_pH_range.shape[1]):
        label = system_names.get(col_idx, f'System {col_idx + 1}')
        # Plot curve
        plt.plot(pH_range_values, Q_pH_range[:, col_idx], label=f'{label} Fit', linewidth=2)
        # Add data point markers
        plt.scatter(pH_range_values, Q_pH_range[:, col_idx], s=60)

    plt.xlabel('pH', fontsize=12)
    plt.ylabel('Q', fontsize=12)
    plt.title(f'Q vs pH at pH 0-14 (Degree {degree} Fit)', fontsize=14)
    plt.legend(fontsize=10, loc='best')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xticks(np.arange(0, 15, 2))
    plt.xlim(0, 14)
    plt.tight_layout()

    # Save range and parameter data
    range_file = os.path.join(plot_dir, 'pH_range_data.txt')
    range_data = np.column_stack([pH_range_values, *Q_pH_range.T])

    # Fix variable name error here (Q_ph_range -> Q_pH_range)
    header = 'pH, ' + ', '.join([f'Q_{system_names.get(i, i + 1)}' for i in range(Q_pH_range.shape[1])])

    np.savetxt(range_file, range_data, header=header, fmt='%.6f')
    print(f"Saved pH range data to: {range_file}")

    param_file = os.path.join(plot_dir, 'fit_params.txt')
    param_header = f'Polynomial parameters (degree {degree}) and R-squared\n'
    param_header += '# Columns: ' + ', '.join(
        [f'coeff_{j} for Poly {i}' for i in range(len(fit_params)) for j in range(degree + 1)[::-1]])
    param_header += ', ' + ', '.join([f'R2_{system_names.get(i, i + 1)}' for i in range(len(fit_params))])
    np.savetxt(param_file, np.column_stack((fit_params, fit_accuracies)), header=param_header, fmt='%.6f')
    print(f"Saved fit parameters to: {param_file}")

    range_plot = os.path.join(plot_dir, 'Q_vs_pH_range.png')
    plt.savefig(range_plot, dpi=300)
    plt.close()
    print(f"Saved Q vs pH range plot: {range_plot}")

    print(f"Q vs pH analysis completed with {len(fit_params)} systems, results saved to {plot_dir}")
def plot_charge_vs_potential(output_dir, file_paths, plot_type="combined",
                             scale_type="SHE", charge_pH=0,DEFAULT_POTENTIAL_POINTS=None):
    """Plot charge vs potential with linear fitting (with SHE/RHE conversion)"""
    # Predefined potential points for charge estimation

    DEFAULT_POTENTIAL_POINTS1=DEFAULT_POTENTIAL_POINTS
    if not file_paths:
        print("Warning: No file paths provided for charge-potential analysis")
        return None

    # Create output directory for charge analysis
    charge_dir = create_output_directory(os.path.join(output_dir, "charge_analysis"))

    # Prepare data structures
    all_data = {}  # Stores main analysis data
    markers = ['o', 's', '^', 'D', 'v']
    colors = plt.cm.viridis(np.linspace(0, 1, len(file_paths)))
    fitting_results = []

    # Define species name mapping
    species_mapping = {i: f"System {i + 1}" for i in range(len(file_paths))}
    if len(file_paths) >= 4:
        species_mapping = {0: 'Slab', 1: '*OOH', 2: '*O', 3: '*OH'}
        if any('o2' in os.path.basename(f).lower() for f in file_paths) and len(file_paths) >= 5:
            species_mapping[4] = '*O₂'

    # Validate scale_type
    valid_scales = ["SHE", "RHE"]
    if scale_type not in valid_scales:
        print(f"Warning: Invalid scale_type '{scale_type}', defaulting to SHE")
        scale_type = "SHE"

    # Validate pH value (only for RHE)
    ph_note = ""  # Global note for figures
    if scale_type == "RHE":
        try:
            charge_pH = float(charge_pH)
            if not (0 <= charge_pH <= 14):
                print(f"Warning: pH={charge_pH} is outside typical range (0-14), using anyway")
            ph_note = f" [pH = {charge_pH:.2f}]"
        except:
            print("Warning: Invalid pH value for RHE conversion, falling back to SHE")
            scale_type = "SHE"
    else:
        # Ignore pH in SHE mode
        charge_pH = 0
        print("Info: pH parameter ignored for SHE scale type")

    print("\n" + "=" * 80)
    print("CHARGE-POTENTIAL ANALYSIS PARAMETERS".center(80))
    print("=" * 80)
    print(f"Voltage Scale: {scale_type}")
    if scale_type == "RHE":
        print(f"pH = {charge_pH:.2f}")
    print("=" * 80)
    print(f"Charge estimation points: {DEFAULT_POTENTIAL_POINTS1} V vs {scale_type}")
    print("=" * 80)

    # SHE to RHE conversion function (corrected formula)
    def U_she_to_U_rhe(U_vs_she, pH):
        """Convert SHE voltage to RHE voltage"""
        return U_vs_she + 0.0592 * pH

    # Process each file
    for idx, file_path in enumerate(file_paths):
        try:
            print(f"Processing file: {file_path}")

            # Read Excel file using read_excel_columns_by_index
            data_columns = read_excel_columns_by_index(file_path)
            if not data_columns or len(data_columns) < 5:
                print(f"File {file_path} has insufficient columns, skipping")
                continue

            # Extract required data columns
            E_fermis = np.array(data_columns[0])
            E_Vacums = np.array(data_columns[1])
            E_potentials = np.array(data_columns[2])

            # Convert data types
            try:
                E_DFTs = np.array([float(x) for x in data_columns[3]])
            except ValueError as e:
                print(f"Error converting DFT Energy data to float in {file_path}: {e}")
                continue

            try:
                Charges = np.array([float(x) for x in data_columns[4]])
            except ValueError as e:
                print(f"Error converting Charge data to float in {file_path}: {e}")
                continue

            # Calculate SHE potential
            _, E_U_vs_SHE = E_free_vs_SHE(E_fermis, E_Vacums, E_potentials, E_DFTs, Charges)

            # Remove invalid points
            valid_mask = np.isfinite(E_U_vs_SHE) & np.isfinite(Charges)
            if np.sum(valid_mask) < 3:
                print(f"Warning: {file_path} insufficient valid data points ({np.sum(valid_mask)}), skipping")
                continue

            U_she = E_U_vs_SHE[valid_mask]
            Q = Charges[valid_mask]

            if len(np.unique(Q)) < 2:
                print(f"Warning: {file_path} has constant charge values, skipping")
                continue

            # Get system name
            system_name = species_mapping.get(idx, os.path.splitext(os.path.basename(file_path))[0])
            safe_name = safe_filename(system_name.replace('*', '_star_'))

            # Convert potential scale if needed
            potential_scale = U_she  # Default SHE scale
            if scale_type == "RHE":
                potential_scale = U_she_to_U_rhe(U_she, charge_pH)
                scale_label = "RHE"
            else:
                scale_label = "SHE"
                ph_note = ""  # Clear pH note for SHE graphs

            # Fit charge vs potential
            slope, intercept = np.polyfit(potential_scale, Q, 1)

            # Calculate R²
            y_pred = slope * potential_scale + intercept
            ss_res = np.sum((Q - y_pred) ** 2)
            ss_tot = np.sum((Q - np.mean(Q)) ** 2)

            if ss_tot < 1e-10:  # Handle case where all points are identical
                r_squared = 0.0
                print(f"Warning for {system_name}: All charge values identical, R² set to 0")
            else:
                r_squared = 1 - (ss_res / ss_tot)

            # Create fit formula
            formula = ""
            additional_info = f"# Linear fit analysis for: {system_name}\n"
            if scale_type == "SHE":
                formula = f"Q = {slope:.4f} × U_SHE + {intercept:.4f}"
                additional_info += f"# Original SHE relation: {formula}\n"
                additional_info += "# REMINDER: Potential points should be vs SHE\n"
            else:  # RHE
                formula = f"Q = {slope:.4f} × U_RHE + {intercept:.4f}"
                # Also show SHE relation for reference
                slope_she_raw, intercept_she_raw = np.polyfit(U_she, Q, 1)
                formula_she_raw = f"Q = {slope_she_raw:.4f} × U_SHE + {intercept_she_raw:.4f}"

                additional_info += f"# Derived relation (RHE): {formula}\n"
                additional_info += f"# Original relation (SHE): {formula_she_raw}\n"
                additional_info += "# REMINDER: Potential points should be vs RHE\n"
                additional_info += f"# pH = {charge_pH} was used for SHE to RHE conversion: U_RHE = U_SHE + 0.0592×pH\n"

            # Store fit results
            fit_result = {
                'system': system_name,
                'scale': scale_label,
                'slope': slope,
                'intercept': intercept,
                'R_squared': r_squared,
                'formula': formula,
                'n_points': len(Q),
                'she_formula': formula_she_raw if scale_type == "RHE" else formula
            }
            fitting_results.append(fit_result)

            # Define interpretation text (physical meaning)
            interpretation = """
# Charge Interpretation (definition):
#
#   Q = (electrons in current system) - (electrons in neutral system)
#
#   • Q > 0: System has a POSITIVE charge (deficit of electrons)
#        ▶ To neutralize: ADD electrons to the system
#        ▶ In VASP: INCREASE NELECT by Q
#
#   • Q < 0: System has a NEGATIVE charge (excess of electrons)
#        ▶ To neutralize: REMOVE electrons from the system
#        ▶ In VASP: DECREASE NELECT by |Q|
#
# How to use at a given potential U:
#   1. Calculate Q = slope × U + intercept
#   2. Adjust NELECT with: NELECT = NELECT_neutral + Q
"""

            # Save raw data and fitting results
            save_file = os.path.join(charge_dir, f"{safe_name}_charge_vs_{scale_label.lower()}.txt")
            image_file = os.path.join(charge_dir, f"charge_vs_{scale_label.lower()}_{safe_name}.png")

            try:
                # Normal write (with Unicode characters)
                with open(save_file, 'w', encoding='utf-8') as f:
                    # Header information
                    f.write(f"# Charge vs {scale_label} Analysis{ph_note}\n")
                    f.write(f"# System: {system_name}\n")
                    f.write(f"# Units: Potential in Volts, Charge in elementary charge units (e)\n\n")

                    # Write additional fit info
                    f.write(additional_info)
                    f.write("\n")

                    # Save data points
                    f.write(f"Potential({scale_label})\tCharge(e)\n")
                    for u, q_val in zip(potential_scale, Q):
                        f.write(f"{u:.6f}\t{q_val:.6f}\n")

                    # Add fitting results
                    f.write("\n# ===== FITTING RESULTS =====\n")
                    f.write(f"# Current scale ({scale_label}): {formula}\n")
                    f.write(f"# Slope (dQ/dU) = {slope:.6f} e/Volt\n")
                    f.write(f"# Intercept (charge at U=0) = {intercept:.6f} e\n")
                    f.write(f"# R-squared (goodness-of-fit) = {r_squared:.6f}\n")
                    f.write(f"# Number of data points = {len(Q)}\n\n")

                    # Add interpretation text
                    f.write(interpretation)

                    # Add RHE conversion info if applicable
                    if scale_label == "RHE":
                        f.write("\n# SHE to RHE conversion details:\n")
                        f.write(f"#   Applied pH: {charge_pH}\n")
                        f.write(f"#   Conversion formula: U_RHE = U_SHE + 0.0592 × pH\n")
                        f.write(
                            "#   Note: pH affects the potential scale but not the intrinsic charge relationship\n\n")

                    # Add charge estimation table
                    f.write(f"\n# Charge Estimation at Key Potentials (vs {scale_label}):\n")
                    f.write("# Potential (V)\tCharge (e)\tNELECT Adjustment\n")

                    for U in DEFAULT_POTENTIAL_POINTS:
                        q_val = slope * U + intercept

                        if q_val > 0:
                            adj_info = f"ADD {q_val:.5f} electrons to neutralize (+{q_val:.5f} e)"
                            vasp_adj = f"Increase NELECT by {q_val:.5f}"
                        elif q_val < 0:
                            adj_info = f"REMOVE {abs(q_val):.5f} electrons to neutralize ({q_val:.5f} e)"
                            vasp_adj = f"Decrease NELECT by {abs(q_val):.5f}"
                        else:
                            adj_info = "Already neutral (no adjustment needed)"
                            vasp_adj = "No change to NELECT"

                        f.write(f"{U:.3f}\t{q_val:.5f}\t{adj_info}\n")
                        f.write(f"#        VASP Operation: {vasp_adj}\n")

            except UnicodeEncodeError:
                # Safe write (ASCII characters only)
                safe_interpretation = interpretation.replace("•", "*").replace("▶", "->")
                with open(save_file, 'w', encoding='utf-8', errors='replace') as f:
                    f.write(f"# Charge vs {scale_label} Analysis{ph_note}\n")
                    f.write(f"# System: {system_name}\n")
                    f.write(f"# Units: Potential in Volts, Charge in elementary charge units (e)\n\n")

                    # Write additional fit info
                    f.write(additional_info)
                    f.write("\n")

                    # Save data points
                    f.write(f"Potential({scale_label})\tCharge(e)\n")
                    for u, q_val in zip(potential_scale, Q):
                        f.write(f"{u:.6f}\t{q_val:.6f}\n")

                    # Add fitting results
                    f.write("\n# ===== FITTING RESULTS =====\n")
                    f.write(f"# Current scale ({scale_label}): {formula}\n")
                    f.write(f"# Slope (dQ/dU) = {slope:.6f} e/Volt\n")
                    f.write(f"# Intercept (charge at U=0) = {intercept:.6f} e\n")
                    f.write(f"# R-squared (goodness-of-fit) = {r_squared:.6f}\n")
                    f.write(f"# Number of data points = {len(Q)}\n\n")

                    # Add interpretation
                    f.write(safe_interpretation)

                    # Add RHE conversion info if applicable
                    if scale_label == "RHE":
                        f.write("\n# SHE to RHE conversion details:\n")
                        f.write(f"#   Applied pH: {charge_pH}\n")
                        f.write(f"#   Conversion formula: U_RHE = U_SHE + 0.0592 × pH\n")
                        f.write(
                            "#   Note: pH affects the potential scale but not the intrinsic charge relationship\n\n")

                    # Add charge estimation table
                    f.write(f"\n# Charge Estimation at Key Potentials (vs {scale_label}):\n")
                    f.write("# Potential (V)\tCharge (e)\tNELECT Adjustment\n")

                    for U in DEFAULT_POTENTIAL_POINTS1:
                        q_val = slope * U + intercept

                        if q_val > 0:
                            adj_info = f"ADD {q_val:.5f} electrons to neutralize (+{q_val:.5f} e)"
                            vasp_adj = f"Increase NELECT by {q_val:.5f}"
                        elif q_val < 0:
                            adj_info = f"REMOVE {abs(q_val):.5f} electrons to neutralize ({q_val:.5f} e)"
                            vasp_adj = f"Decrease NELECT by {abs(q_val):.5f}"
                        else:
                            adj_info = "Already neutral (no adjustment needed)"
                            vasp_adj = "No change to NELECT"

                        f.write(f"{U:.3f}\t{q_val:.5f}\t{adj_info}\n")
                        f.write(f"#        VASP Operation: {vasp_adj}\n")

            print(f"Saved {scale_label} data for {system_name} to: {save_file}")
            all_data[system_name] = (potential_scale, Q, slope, intercept, scale_label)

            # Create individual plot if requested
            if plot_type == "individual":
                fig, ax = plt.subplots(figsize=(9, 7))

                # Scatter plot
                ax.scatter(potential_scale, Q, marker=markers[idx % len(markers)],
                           color=colors[idx], s=100, alpha=0.8,
                           label=f'{system_name} Data')

                # Fit line
                fit_x = np.linspace(min(potential_scale), max(potential_scale), 100)
                fit_y = slope * fit_x + intercept
                ax.plot(fit_x, fit_y, 'r-', linewidth=2.5,
                        label=f'Fit: {formula}')

                # Zero references
                ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
                ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)

                # Set labels and title
                ax.set_xlabel(f'Potential vs {scale_label} (V)', fontsize=14)
                ax.set_ylabel('Charge (e)', fontsize=14)
                title = f'Charge vs Potential: {system_name}\n{formula} (R² = {r_squared:.4f})' + ph_note
                ax.set_title(title, fontsize=16)
                ax.grid(True, linestyle='--', alpha=0.5)

                # Add annotation box
                annotation_text = f"Slope: {slope:.4f} e/V\nR²: {r_squared:.4f}"
                ax.annotate(annotation_text, xy=(0.05, 0.95), xycoords='axes fraction',
                            fontsize=12, bbox=dict(boxstyle="round,pad=0.3", fc='white', alpha=0.8))

                ax.legend(loc='best', fontsize=12)
                plt.tight_layout()

                # Save image
                plt.savefig(image_file, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"Saved {scale_label} plot for {system_name} to: {image_file}")

        except Exception as e:
            print(f"Error processing file {file_path}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue

    # Combined plot
    if plot_type == "combined" and all_data:
        fig, ax = plt.subplots(figsize=(12, 9))
        plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
        plt.axvline(x=0, color='k', linestyle='--', alpha=0.3)

        for idx, (system_name, data) in enumerate(all_data.items()):
            potential_scale, Q, slope, intercept, scale_label = data
            color = colors[idx % len(colors)]

            # Scatter plot
            ax.scatter(potential_scale, Q, marker=markers[idx % len(markers)],
                       color=color, s=120, alpha=0.8,
                       label=f'{system_name} Data')

            # Fit line
            fit_x = np.linspace(min(potential_scale) - 0.1, max(potential_scale) + 0.1, 100)
            fit_y = slope * fit_x + intercept
            ax.plot(fit_x, fit_y, color=color, linestyle='-', linewidth=2.5,
                    label=f'{system_name}: {fitting_results[idx]["formula"]}')

        # Set labels and title
        ax.set_xlabel(f'Potential vs {scale_label} (V)', fontsize=14)
        ax.set_ylabel('Charge (e)', fontsize=14)
        title = f'Charge vs {scale_label} Potential' + ph_note
        ax.set_title(title, fontsize=16)
        ax.grid(True, linestyle='--', alpha=0.5)

        # Adjust legend
        ax.legend(loc='best', fontsize=10, ncol=2)
        plt.tight_layout()

        # Save combined plot
        image_file = os.path.join(charge_dir, f'charge_vs_{scale_label.lower()}_combined.png')
        plt.savefig(image_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved combined {scale_label} plot to: {image_file}")

    # Save all fitting results
    if fitting_results:
        fit_file = os.path.join(charge_dir, f"charge_fitting_results_{scale_label.lower()}.csv")

        try:
            with open(fit_file, 'w', encoding='utf-8') as f:
                f.write(
                    "System,Potential Scale,dQ/dU (e/Volt),Intercept (e),R_squared,Formula,She_Formula,N_points,Electron Adjustment Rule\n")
                for res in fitting_results:
                    she_formula = res.get('she_formula', res['formula'])
                    slope = res['slope']
                    intercept = res['intercept']
                    rule = "NELECT = NELECT_neutral + Q"

                    f.write(
                        f"{res['system']},{res['scale']},{slope:.6f},{intercept:.6f},{res['R_squared']:.6f},\"{res['formula']}\",\"{she_formula}\",{res['n_points']},\"{rule}\"\n")

        except UnicodeEncodeError:
            # Safe write for CSV
            with open(fit_file, 'w', encoding='utf-8', errors='replace') as f:
                f.write(
                    "System,Potential Scale,dQ/dU (e/Volt),Intercept (e),R_squared,Formula,She_Formula,N_points,Electron Adjustment Rule\n")
                for res in fitting_results:
                    she_formula = res.get('she_formula', res['formula'])
                    slope = res['slope']
                    intercept = res['intercept']
                    rule = "NELECT = NELECT_neutral + Q"

                    f.write(
                        f"{res['system']},{res['scale']},{slope:.6f},{intercept:.6f},{res['R_squared']:.6f},\"{res['formula']}\",\"{she_formula}\",{res['n_points']},\"{rule}\"\n")

        print(f"Saved all fitting results to: {fit_file}")

        # Print summary
        print("\n\n" + "=" * 90)
        print(" CHARGE-POTENTIAL RELATIONSHIP SUMMARY ".center(90))
        print("=" * 90)
        print(f"{'System':<15}{'Fit Equation':<40}{'R²':<8}{'Data Points':>10}")
        print("-" * 90)
        for res in fitting_results:
            scale_note = f" (vs {res['scale']})"
            print(f"{res['system']:<15}{res['formula']:<40}{res['R_squared']:.4f}{'':<5}{res['n_points']:>5}")
            if res.get('she_formula') and res['scale'] == "RHE":
                print(f"{' ':>15}SHE Form: {res['she_formula']}")
        print("=" * 90)

        # Physical interpretation and VASP guide
        print("\nPHYSICAL INTERPRETATION & VASP APPLICATION GUIDE:")
        print(" ■ Charge definition: Q = (electrons in system) - (electrons in neutral system)")
        print(f" ■ Potential scale: {scale_type} (as selected) - ALWAYS use this scale when applying formulas")
        if scale_type == "RHE":
            print(f"    → pH = {charge_pH:.2f} was used for conversion: U_RHE = U_SHE + 0.0592 × pH")

        print("\n ■ Charge sign conventions:")
        print("    → Q > 0: POSITIVE charge (electron deficit)")
        print("        ▶ In VASP: INCREASE NELECT by Q (+Q e)")
        print("        ▶ Physics: System needs ADDITIONAL electrons to become neutral")
        print("    → Q < 0: NEGATIVE charge (electron excess)")
        print("        ▶ In VASP: DECREASE NELECT by |Q| (-|Q| e)")
        print("        ▶ Physics: System needs to REMOVE electrons to become neutral")

        print("\n ■ How to use at a specific potential:")
        print(f"    1. Choose potential U (vs {scale_type}) from: {DEFAULT_POTENTIAL_POINTS1} or your desired value")
        print("    2. Calculate Q = slope × U + intercept (using the fit equation)")
        print("    3. Adjust NELECT using: NELECT = NELECT_neutral + Q")

        print("\n ■ Important notes:")
        if scale_type == "RHE":
            print("    ・ For RHE conversions, maintain the same pH in simulations")
            print("    ・ Your input POTIMPA file MUST contain SHE potentials (they will be converted to RHE)")
        else:
            print("    ・ Your input POTIMPA file MUST contain SHE potentials")

        print("    ・ VASP NELECT must be specified as: NELECT = NELECT_neutral + Q")
        print("    ・ Accurately determine NELECT_neutral for your neutral system")
        print("    ・ For adspecies calculations, consider coupled proton-electron transfers")

        print("=" * 90)
        print(f"\n[SUCCESS] Charge analysis completed against {scale_type}. Results saved to: {charge_dir}")

    return fitting_results
def plot_ads_energy_vs_U(output_dir, ads_systems, fit_params_list, voltage_scale="SHE", pH=0, U_range=None):
    """Plot adsorption energy vs U (SHE or RHE) with customizable U range"""
    plot_dir = create_output_directory(os.path.join(output_dir, "ads_energy_plots"))
    plt.figure(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(ads_systems)))

    if U_range is None:
        U_range = np.linspace(-2, 1, 100)  # Default range if not provided

    for i, system in enumerate(ads_systems):
        C, U_pzc, E0, _ = fit_params_list[i]

        if voltage_scale == "RHE":
            U_range_she = U_rhe_to_U_she(U_range, pH)
        else:
            U_range_she = U_range

        E_ads = [calculate_ads_energy(u, C, U_pzc, E0) for u in U_range_she]

        plt.plot(U_range, E_ads, color=colors[i], label=system, linewidth=2)

    plt.xlabel(f'Potential vs {voltage_scale} (V)', fontsize=12)
    plt.ylabel('Adsorption Energy ΔE (eV)', fontsize=12)
    plt.title('Adsorption Energy vs Potential', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f"ads_energy_vs_{voltage_scale}_pH{pH}.png"), dpi=300)
    plt.close()

def plot_E_free_vs_U_she(output_dir, data_source, data, file_paths, R2_list=None):
    """Plot free energy vs SHE potential fitting curve and output raw data (with R² values)"""
    plot_dir = create_output_directory(os.path.join(output_dir, "free_energy_fits"))
    raw_data_dir = create_output_directory(os.path.join(plot_dir, "raw_fit_data"))

    # Initialize plot
    plt.figure(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, 10))

    # Store all fitting data
    all_fit_data = {}

    # Check if R2 list exists
    if R2_list is None:
        R2_list = []

    # ======== Process literature parameters ========
    if data_source == "literature_params":
        species_order = ["slab", "*OOH", "*O", "*OH"]

        for i, species in enumerate(species_order):
            if species not in data:
                print(f"Warning: Species {species} not found in literature parameters")
                continue

            try:
                params = [
                    data[species]["C"],
                    data[species]["U_pzc"],
                    data[species]["E0"]
                ]
                C, U_pzc, E0 = params

                # Create voltage range
                min_U = min(-3, U_pzc - 1.5)
                max_U = max(3, U_pzc + 1.5)
                x_range = np.linspace(min_U, max_U, 100)
                y_fit = quadratic_function_vs_she(x_range, C, U_pzc, E0)

                # Use current color
                current_color = colors[i % len(colors)]
                plt.plot(x_range, y_fit, color=current_color, linewidth=2, label=species)

                # Create safe name
                safe_name = species.replace('*', '').replace(' ', '_')
                all_fit_data[safe_name] = {
                    "fitted": np.column_stack((x_range, y_fit)),
                    "params": params,
                    "R2": None  # Literature parameters don't have R²
                }

                # Create subdirectory to save data
                species_dir = os.path.join(raw_data_dir, safe_name)
                os.makedirs(species_dir, exist_ok=True)

                # Save fitted data
                fitted_out = os.path.join(species_dir, "fitted_data.txt")
                np.savetxt(fitted_out, all_fit_data[safe_name]["fitted"],
                           header="U_range\tFit_energy", fmt="%.8f")

                # Save fitting parameters with utf-8 encoding
                param_out = os.path.join(species_dir, "fit_params.txt")
                with open(param_out, "w", encoding='utf-8') as f:
                    f.write(f"System: {species}\n")
                    f.write(f"Coefficient (C): {C:.8f}\n")
                    f.write(f"U_pzc: {U_pzc:.8f}\n")
                    f.write(f"E0: {E0:.8f}\n")
                    f.write("R²: Not available (pre-fitted literature parameters)\n")

                print(f"Saved {species} fitting data to {species_dir}")

            except Exception as e:
                print(f"Error processing literature parameter {species}: {str(e)}")

    # Process raw literature data
    # Process literature data
    elif data_source == "literature":
        systems = list(data.keys())

        for i, system in enumerate(systems):
            try:
                U = np.array(data[system]["U"])
                E = np.array(data[system]["E"])
                if len(U) < 3 or len(E) < 3:
                    print(f"Warning: {system} has insufficient data points ({len(U)}), skipping")
                    continue

                # Remove invalid points
                mask = np.isfinite(U) & np.isfinite(E)
                U = U[mask]
                E = E[mask]

                if len(U) < 3:
                    print(f"Warning: {system} has insufficient valid data points ({len(U)}), skipping")
                    continue

                # Dynamic initial guess
                C_init = 1.0  # Typical capacitance-like coefficient
                U_pzc_init = np.mean(U) if len(U) > 0 else 0.0  # Center of potential range
                E0_init = np.mean(E) if len(E) > 0 else 0.0  # Center of energy range

                # Parameter bounds: C > 0, U_pzc within reasonable voltage, E0 within reasonable energy
                bounds = ([0, -5, -1000], [10, 5, 1000])  # Adjust based on expected ranges

                # Perform fitting with bounds and increased maxfev
                params, _ = curve_fit(
                    quadratic_function_vs_she,
                    U,
                    E,
                    p0=[C_init, U_pzc_init, E0_init],
                    bounds=bounds,
                    maxfev=10000,  # Increase max function evaluations
                    method='trf'  # Trust Region Reflective for robustness
                )
                C, U_pzc, E0 = params

                # Calculate R²
                E_pred = quadratic_function_vs_she(U, C, U_pzc, E0)
                ss_res = np.sum((E - E_pred) ** 2)
                ss_tot = np.sum((E - np.mean(E)) ** 2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else float('nan')

                # Create fitted curve
                x_range = np.linspace(min(U) - 0.5, max(U) + 0.5, 100)
                y_fit = quadratic_function_vs_she(x_range, C, U_pzc, E0)

                # Get current color
                current_color = colors[i % len(colors)]
                plt.plot(x_range, y_fit, color=current_color, linewidth=2)
                plt.scatter(U, E, color=current_color, s=60, marker='o', label=system)

                # Create safe name
                safe_name = safe_filename(system.replace('*', '_star_'))
                all_fit_data[safe_name] = {
                    "original": np.column_stack((U, E)),
                    "fitted": np.column_stack((x_range, y_fit)),
                    "params": params,
                    "R2": r_squared
                }

                # Create subdirectory
                system_dir = os.path.join(raw_data_dir, safe_name)
                os.makedirs(system_dir, exist_ok=True)

                # Save raw data (English, utf8 encoding)
                original_out = os.path.join(system_dir, "original_data.txt")
                with open(original_out, 'w', encoding='utf-8') as f:
                    f.write("# U_calc\tE_calc\n")
                    for u_val, e_val in zip(U, E):
                        f.write(f"{u_val}\t{e_val}\n")

                # Save fitted data (English, utf8 encoding)
                fitted_out = os.path.join(system_dir, "fitted_data.txt")
                with open(fitted_out, 'w', encoding='utf-8') as f:
                    f.write("# U_range\tFit_energy\n")
                    for x, y in zip(x_range, y_fit):
                        f.write(f"{x}\t{y}\n")

                # Save fitting parameters (English, utf8 encoding)
                param_out = os.path.join(system_dir, "fit_params.txt")
                with open(param_out, 'w', encoding='utf-8') as f:
                    f.write(f"System: {system}\n")
                    f.write(f"Coefficient (C): {C:.8f}\n")
                    f.write(f"U_pzc: {U_pzc:.8f}\n")
                    f.write(f"E0: {E0:.8f}\n")
                    f.write(f"Coefficient of Determination (R squared): {r_squared:.8f}\n")

                print(f"{system} fitted! C={C:.4f}, U_pzc={U_pzc:.4f}, E0={E0:.4f}, R2={r_squared:.4f}")
                print(f"Data saved to: {system_dir}")

            except Exception as e:
                print(f"Error fitting {system}: {str(e)}")

    # Process file data
    elif data_source == "file" and file_paths:
        for idx, file_path in enumerate(file_paths):
            try:
                # Read Excel column data
                data_columns = read_excel_columns_by_index(file_path)
                if not data_columns or len(data_columns) < 5:
                    print(f"File {file_path} format error, skipping")
                    continue

                # Calculate free energy
                E_fermis = np.array(data_columns[0])
                E_Vacums = np.array(data_columns[1])
                E_potentials = np.array(data_columns[2])
                E_DFTs = np.array(data_columns[3])
                Charges = np.array(data_columns[4])

                E_free_U, E_U_vs_SHE = E_free_vs_SHE(
                    E_fermis, E_Vacums, E_potentials, E_DFTs, Charges)

                # Remove invalid points
                mask = np.isfinite(E_U_vs_SHE) & np.isfinite(E_free_U)
                E_U_vs_SHE = E_U_vs_SHE[mask]
                E_free_U = E_free_U[mask]

                if len(E_U_vs_SHE) < 3:
                    print(f"{file_path} insufficient valid points ({len(E_U_vs_SHE)}), skipping")
                    continue

                # Dynamic initial guess
                C_init = 1.0  # Typical capacitance-like coefficient
                U_pzc_init = np.mean(E_U_vs_SHE) if len(E_U_vs_SHE) > 0 else 0.0  # Center of potential range
                E0_init = np.mean(E_free_U) if len(E_free_U) > 0 else 0.0  # Center of energy range

                # Parameter bounds: C > 0, U_pzc within reasonable voltage, E0 within reasonable energy
                bounds = ([0, -5, -1000], [10, 5, 1000])  # Adjust based on expected ranges

                # Perform fitting with bounds and increased maxfev
                params, _ = curve_fit(
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

                # Create fitted curve
                x_range = np.linspace(min(E_U_vs_SHE) - 0.5, max(E_U_vs_SHE) + 0.5, 100)
                y_fit = quadratic_function_vs_she(x_range, C, U_pzc, E0)

                # Get current color
                current_color = colors[idx % len(colors)]
                plt.plot(x_range, y_fit, color=current_color, linewidth=2)

                # Get filename for label
                file_name = os.path.splitext(os.path.basename(file_path))[0]
                plt.scatter(E_U_vs_SHE, E_free_U, color=current_color, s=60, marker='o', label=file_name)

                safe_name = safe_filename(file_name)
                all_fit_data[safe_name] = {
                    "original": np.column_stack((E_U_vs_SHE, E_free_U)),
                    "fitted": np.column_stack((x_range, y_fit)),
                    "params": params,
                    "R2": r_squared
                }

                # Create subdirectory
                file_dir = os.path.join(raw_data_dir, safe_name)
                os.makedirs(file_dir, exist_ok=True)

                # Save raw data (English, utf8 encoding)
                original_out = os.path.join(file_dir, "original_data.txt")
                with open(original_out, 'w', encoding='utf-8') as f:
                    f.write("# U_calc\tE_calc\n")
                    for u_val, e_val in zip(E_U_vs_SHE, E_free_U):
                        f.write(f"{u_val}\t{e_val}\n")

                # Save fitted data (English, utf8 encoding)
                fitted_out = os.path.join(file_dir, "fitted_data.txt")
                with open(fitted_out, 'w', encoding='utf-8') as f:
                    f.write("# U_range\tFit_energy\n")
                    for x, y in zip(x_range, y_fit):
                        f.write(f"{x}\t{y}\n")

                # Save fitting parameters (English, utf8 encoding)
                param_out = os.path.join(file_dir, "fit_params.txt")
                with open(param_out, 'w', encoding='utf-8') as f:
                    f.write(f"System: {file_name}\n")
                    f.write(f"Coefficient (C): {C:.8f}\n")
                    f.write(f"U_pzc: {U_pzc:.8f}\n")
                    f.write(f"E0: {E0:.8f}\n")
                    f.write(f"Coefficient of Determination (R squared): {r_squared:.8f}\n")

                print(f"{file_name}: C={C:.4f}, U_pzc={U_pzc:.4f}, E0={E0:.4f}, R2={r_squared:.4f}")
                print(f"Data saved to: {file_dir}")

            except Exception as e:
                print(f"Error processing file {file_path}: {str(e)}")

    # Set plot properties (English)
    plt.xlabel('Potential vs SHE (V)', fontsize=12)
    plt.ylabel('Free Energy (eV)', fontsize=12)
    plt.title('Free Energy vs SHE Potential', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()

    # Save plot (English filename)
    plt.savefig(os.path.join(plot_dir, "free_energy_vs_she.png"), dpi=300)
    plt.close()

    # ===== Save all data to a combined file (English content) =====
    if all_fit_data:
        all_fit_file = os.path.join(raw_data_dir, "all_fit_data.csv")
        with open(all_fit_file, 'w', encoding='utf-8') as f:
            f.write("System,U,Free_Energy\n")
            for safe_name, data_dict in all_fit_data.items():
                if "original" in data_dict:
                    data_to_write = data_dict["original"]
                else:
                    data_to_write = data_dict["fitted"]

                # Iterate all data points
                for row in data_to_write:
                    f.write(f"{safe_name},{row[0]},{row[1]}\n")

    # Create parameter summary file (English content)
    if all_fit_data:
        summary_file = os.path.join(raw_data_dir, "fit_parameters_summary.csv")
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("System,C,U_pzc,E0,R_squared\n")
            for safe_name, data_dict in all_fit_data.items():
                params = data_dict["params"]
                r2 = data_dict.get("R2", "N/A")
                f.write(f"{safe_name},{params[0]},{params[1]},{params[2]},{r2}\n")

        print(f"Summary saved to: {summary_file}")

    print(f"Analysis completed! All results saved to: {raw_data_dir}")

def plot_dG_diagrams(output_dir, results, pH, reaction_type, voltage_scale):
    voltage_range = results["voltage"]
    intermediates = results["intermediates"]
    steps = results["steps"]

    create_output_directory(output_dir)

    plt.figure(figsize=(10, 6))
    plt.plot(voltage_range, intermediates["delta_GOOH"], label='*OOH', linewidth=2)
    plt.plot(voltage_range, intermediates["delta_GO"], label='*O', linewidth=2)
    plt.plot(voltage_range, intermediates["delta_GOH"], label='*OH', linewidth=2)

    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    plt.xlabel(f'Potential vs {voltage_scale} (V)', fontsize=12)
    plt.ylabel('ΔG (eV)', fontsize=12)
    plt.title(f'Intermediates Energy (pH={pH}, {reaction_type})', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'intermediates_pH{pH}_{reaction_type}_VS_{voltage_scale}.png'), dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(voltage_range, steps["d_G1"], label='Step 1', linewidth=2)
    plt.plot(voltage_range, steps["d_G2"], label='Step 2', linewidth=2)
    plt.plot(voltage_range, steps["d_G3"], label='Step 3', linewidth=2)
    plt.plot(voltage_range, steps["d_G4"], label='Step 4', linewidth=2)

    plt.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.7)
    plt.xlabel(f'Potential vs {voltage_scale} (V)', fontsize=12)
    plt.ylabel('ΔG (eV)', fontsize=12)
    plt.title(f'Step Energy (pH={pH}, {reaction_type})', fontsize=14)
    plt.legend(fontsize=9)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'steps_pH{pH}_{reaction_type}_VS_{voltage_scale}.png'), dpi=300)
    plt.close()

def plot_heatmap(X, Y, Z, title, xlabel, ylabel, zlabel, filename):
    try:
        plt.figure(figsize=(10, 8))

        if not np.any(np.isfinite(Z)):
            print(f"Warning: No valid data for {title}, skipping plot")
            return

        z_min, z_max = np.nanmin(Z), np.nanmax(Z)

        colors = [
            (0.0, '#313695'),
            (0.25, '#4575b4'),
            (0.4, '#74add1'),
            (0.55, '#a8ddb5'),
            (0.7, '#fee08b'),
            (0.85, '#fdae61'),
            (1.0, '#f46d43')
        ]
        cmap = LinearSegmentedColormap.from_list('diverging', colors)

        X_grid, Y_grid = np.meshgrid(X, Y)

        pc = plt.pcolormesh(X_grid, Y_grid, Z, cmap=cmap, shading='auto', vmin=z_min, vmax=z_max)

        valid_z = Z[np.isfinite(Z)]
        if len(valid_z) > 0:
            contour_levels = np.arange(z_min, z_max, (z_max - z_min) / 15)
            plt.contour(X_grid, Y_grid, Z, levels=contour_levels, colors='k', linewidths=0.5, alpha=0.7)

        cbar = plt.colorbar(pc, pad=0.1, aspect=30, format='%.2f')
        cbar.set_label(zlabel, rotation=90, labelpad=15, fontsize=12)

        plt.xlabel(xlabel, fontsize=12)
        plt.ylabel(ylabel, fontsize=12)
        plt.title(title, fontsize=14, pad=12)

        if len(X) > 1 and max(X) - min(X) >= 1:
            plt.xticks(np.arange(min(X), max(X) + 0.5, 1), fontsize=10)
        else:
            plt.xticks(fontsize=10)

        if len(Y) > 1 and max(Y) - min(Y) >= 0.1:
            plt.yticks(np.arange(min(Y), max(Y) + 0.1, 0.2), fontsize=10)
        else:
            plt.yticks(fontsize=10)

        plt.xlim(min(X), max(X))
        plt.ylim(min(Y), max(Y))

        plt.grid(True, linestyle='--', alpha=0.4)

        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Heatmap saved: {filename}")

    except Exception as e:
        print(f"Error plotting heatmap '{title}': {str(e)}")