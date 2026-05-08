# calculations.py
import numpy as np
from cp_analyzer.config_default import G_O2, G_H2, G_H2O, delta_G_water
from cp_analyzer.utils import U_rhe_to_U_she, quadratic_function_vs_she, U_she_to_U_rhe

import os
import json
import pandas as pd
from cp_analyzer.data_processing import create_output_directory
from cp_analyzer.plotting import plot_heatmap
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, brentq
def calculate_intermediates(voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2=True):
    intermediates = []
    for i in range(len(C_list)):
        E_vs = quadratic_function_vs_she(voltage_she, C_list[i], U_pzc_list[i], E0_list[i])
        intermediates.append(E_vs + E_zpe_st[i])

    if include_O2 and len(intermediates) >= 5:
        G_sub, G_oo, G_ooh, G_o, G_oh = intermediates[:5]
    elif len(intermediates) >= 4:
        G_sub, G_ooh, G_o, G_oh = intermediates[:4]
        G_oo = G_sub + G_O2
    else:
        raise ValueError("Insufficient intermediates data")

    delta_GOOH = G_ooh - G_sub - (2 * G_H2O - 1.5 * G_H2)
    delta_GO = G_o - G_sub - (G_H2O - G_H2)
    delta_GOH = G_oh - G_sub - (G_H2O - 0.5 * G_H2)

    return delta_GOOH, delta_GO, delta_GOH

def calculate_steps_from_intermediates(delta_GOOH, delta_GO, delta_GOH, voltage_rhe, reaction_type="ORR"):
    if reaction_type == "ORR":
        d_G1 = delta_GOOH - delta_G_water + voltage_rhe
        d_G2 = delta_GO - delta_GOOH + voltage_rhe
        d_G3 = delta_GOH - delta_GO + voltage_rhe
        d_G4 = -delta_GOH + voltage_rhe
        return d_G1, d_G2, d_G3, d_G4
    else:  # OER
        d_G1 = delta_GOH - voltage_rhe
        d_G2 = (delta_GO - delta_GOH) - voltage_rhe
        d_G3 = (delta_GOOH - delta_GO) - voltage_rhe
        d_G4 = (delta_G_water - delta_GOOH) - voltage_rhe
        return d_G1, d_G2, d_G3, d_G4


def calculate_for_voltage_range(voltage_range, pH, C_list, U_pzc_list, E0_list, E_zpe_st,
                                reaction_type, include_O2=True, voltage_scale="RHE"):
    """
    Calculate ORR/OER free energy for a range of voltages at given pH.

    Parameters:
    - voltage_range: array of voltages (unit depends on voltage_scale)
    - voltage_scale: "RHE" or "SHE" — determines how input voltages are interpreted
    """
    if len(voltage_range) == 0:
        return {
            "voltage": np.array([]),
            "intermediates": {
                "delta_GOOH": np.array([]),
                "delta_GO": np.array([]),
                "delta_GOH": np.array([])
            },
            "steps": {
                "d_G1": np.array([]),
                "d_G2": np.array([]),
                "d_G3": np.array([]),
                "d_G4": np.array([])
            }
        }

    delta_GOOH_list = []
    delta_GO_list = []
    delta_GOH_list = []
    d_G1_list = []
    d_G2_list = []
    d_G3_list = []
    d_G4_list = []

    print(f"Calculating {reaction_type} at pH={pH}, voltage_scale={voltage_scale}, "
          f"voltage range: {min(voltage_range):.3f} ~ {max(voltage_range):.3f} V")

    for voltage_input in voltage_range:
        try:
            # ==================== 电压转换逻辑（核心修改部分） ====================
            if voltage_scale == "RHE":
                # 输入是 RHE → 转为 SHE 用于中间体能量计算（拟合基于 SHE）
                voltage_she = U_rhe_to_U_she(voltage_input, pH)
                voltage_for_steps = voltage_input  # 步骤计算使用 RHE
            else:  # SHE
                # 输入是 SHE → 直接使用
                voltage_she = voltage_input
                voltage_for_steps = voltage_input  # 步骤计算也使用 SHE

            # 计算中间体吸附自由能（必须用 SHE 电压）
            delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
                voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2
            )

            # 计算各步反应自由能
            dG1, dG2, dG3, dG4 = calculate_steps_from_intermediates(
                delta_GOOH, delta_GO, delta_GOH,
                voltage_for_steps, reaction_type
            )

            # 保存结果
            delta_GOOH_list.append(delta_GOOH)
            delta_GO_list.append(delta_GO)
            delta_GOH_list.append(delta_GOH)
            d_G1_list.append(dG1)
            d_G2_list.append(dG2)
            d_G3_list.append(dG3)
            d_G4_list.append(dG4)

        except Exception as e:
            print(f"Warning: Error at voltage={voltage_input:.3f} V ({voltage_scale}), pH={pH}: {str(e)}")
            delta_GOOH_list.append(float('nan'))
            delta_GO_list.append(float('nan'))
            delta_GOH_list.append(float('nan'))
            d_G1_list.append(float('nan'))
            d_G2_list.append(float('nan'))
            d_G3_list.append(float('nan'))
            d_G4_list.append(float('nan'))

    results = {
        "voltage": np.array(voltage_range),
        "intermediates": {
            "delta_GOOH": np.array(delta_GOOH_list),
            "delta_GO": np.array(delta_GO_list),
            "delta_GOH": np.array(delta_GOH_list)
        },
        "steps": {
            "d_G1": np.array(d_G1_list),
            "d_G2": np.array(d_G2_list),
            "d_G3": np.array(d_G3_list),
            "d_G4": np.array(d_G4_list)
        }
    }

    # 验证：ORR 在 U=0 RHE 时四步之和应接近 -4.92 eV
    if reaction_type == "ORR" and voltage_scale == "RHE":
        idx_zero = np.argmin(np.abs(voltage_range - 0.0))
        total_dG = (results["steps"]["d_G1"][idx_zero] +
                    results["steps"]["d_G2"][idx_zero] +
                    results["steps"]["d_G3"][idx_zero] +
                    results["steps"]["d_G4"][idx_zero])
        print(f"  Verification at U≈0 RHE: sum(dG1~dG4) = {total_dG:.4f} eV  (should be ≈ -4.92)")

    return results
def calculate_dG_grid(pH_range, voltage_range, C_list, U_pzc_list, E0_list, E_zpe_st, reaction_type="ORR", include_O2=True, voltage_scale="RHE"):
    n_pH = len(pH_range)
    n_voltage = len(voltage_range)

    print(f"Calculating {reaction_type} grid: {n_pH}x{n_voltage} points")
    print(f"pH range: {min(pH_range)}-{max(pH_range)}")
    print(f"Voltage range ({voltage_scale}): {min(voltage_range)}-{max(voltage_range)}V")

    grid_data = {
        "reaction": reaction_type,
        "voltage_scale": voltage_scale,
        "intermediates": {
            "delta_GOOH": np.zeros((n_pH, n_voltage)),
            "delta_GO": np.zeros((n_pH, n_voltage)),
            "delta_GOH": np.zeros((n_pH, n_voltage))
        },
        "steps": {
            "d_G1": np.zeros((n_pH, n_voltage)),
            "d_G2": np.zeros((n_pH, n_voltage)),
            "d_G3": np.zeros((n_pH, n_voltage)),
            "d_G4": np.zeros((n_pH, n_voltage))
        },
        "pH": np.array(pH_range),
        "voltage": np.array(voltage_range)
    }

    for i, pH_val in enumerate(pH_range):
        if (i + 1) % 5 == 0 or (i + 1) == n_pH:
            print(f"  Progress: {i + 1}/{n_pH} pH points...")

        for j, voltage_val in enumerate(voltage_range):
            if voltage_scale == "RHE":
                voltage_she = U_rhe_to_U_she(voltage_val, pH_val)
            else:
                voltage_she = voltage_val

            try:
                delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
                    voltage_she, pH_val, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2
                )

                d_G1, d_G2, d_G3, d_G4 = calculate_steps_from_intermediates(
                    delta_GOOH, delta_GO, delta_GOH, voltage_val, reaction_type
                )

                grid_data["intermediates"]["delta_GOOH"][i, j] = delta_GOOH
                grid_data["intermediates"]["delta_GO"][i, j] = delta_GO
                grid_data["intermediates"]["delta_GOH"][i, j] = delta_GOH
                grid_data["steps"]["d_G1"][i, j] = d_G1
                grid_data["steps"]["d_G2"][i, j] = d_G2
                grid_data["steps"]["d_G3"][i, j] = d_G3
                grid_data["steps"]["d_G4"][i, j] = d_G4
            except Exception as e:
                print(f"Error calculating at pH={pH_val}, U={voltage_val}: {e}")

    return grid_data

def save_dG_grid(grid_data, output_dir):
    save_dir = create_output_directory(os.path.join(output_dir, "grid_data"))
    prefix = grid_data["reaction"]

    intermediates = grid_data["intermediates"]
    for key, data in intermediates.items():
        filename = f"{prefix}_intermediate_{key}.csv"
        df = pd.DataFrame(data, index=grid_data["pH"], columns=grid_data["voltage"])
        df.index.name = 'pH'
        df.columns.name = f'Voltage ({grid_data["voltage_scale"]})'
        df.to_csv(os.path.join(save_dir, filename))

    steps = grid_data["steps"]
    for key, data in steps.items():
        filename = f"{prefix}_step_{key}.csv"
        df = pd.DataFrame(data, index=grid_data["pH"], columns=grid_data["voltage"])
        df.index.name = 'pH'
        df.columns.name = f'Voltage ({grid_data["voltage_scale"]})'
        df.to_csv(os.path.join(save_dir, filename))

    print(f"Grid data saved to: {save_dir}")

def plot_dG_heatmaps(grid_data, output_dir):
    plot_dir = create_output_directory(os.path.join(output_dir, "heatmaps", grid_data["reaction"]))
    prefix = grid_data["reaction"]

    pH_range = grid_data["pH"]
    voltage_range = grid_data["voltage"]

    X, Y = np.meshgrid(pH_range, voltage_range)

    intermediates = grid_data["intermediates"]
    for key, data in intermediates.items():
        if np.isnan(data).all():
            print(f"Warning: No valid data for {prefix} {key}, skipping plot")
            continue

        data_transposed = data.T

        plot_heatmap(
            pH_range, voltage_range, data_transposed,
            f'{prefix} Intermediate: {key}',
            'pH',
            f'Potential (V vs {grid_data["voltage_scale"]})',
            f'ΔG_{key} (eV)',
            os.path.join(plot_dir, f'heatmap_{prefix}_{key}.png')
        )

    steps = grid_data["steps"]
    for key, data in steps.items():
        if np.isnan(data).all():
            print(f"Warning: No valid data for {prefix} {key}, skipping plot")
            continue

        data_transposed = data.T

        plot_heatmap(
            pH_range, voltage_range, data_transposed,
            f'{prefix} Step: {key}',
            'pH',
            f'Potential (V vs {grid_data["voltage_scale"]})',
            f'ΔG_{key} (eV)',
            os.path.join(plot_dir, f'heatmap_{prefix}_{key}.png')
        )

    print(f"{prefix} heatmap analysis completed")


def save_microkinetics_parameters(output_dir, pH_values, voltage_values, C_list, U_pzc_list, E0_list, E_zpe_st,
                                  include_O2=True, voltage_scale="RHE"):
    save_dir = create_output_directory(os.path.join(output_dir, "microkinetics"))
    prefix = f"microkinetics_{voltage_scale}"

    print(f"Generating MICROKINETICS parameters (pH points: {len(pH_values)}, potential points: {len(voltage_values)})")

    orr_data = []
    oer_data = []

    for pH in pH_values:
        for voltage in voltage_values:
            try:
                if voltage_scale == "RHE":
                    voltage_she = U_rhe_to_U_she(voltage, pH)
                else:
                    voltage_she = voltage

                delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
                    voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2
                )

                d_orr = calculate_steps_from_intermediates(delta_GOOH, delta_GO, delta_GOH, voltage,"ORR")

                d_oer = calculate_steps_from_intermediates(delta_GOOH, delta_GO, delta_GOH, voltage,"OER")

                orr_entry = f"{pH:.4f}, {voltage:.4f}, {d_orr[0]:.6f}, {d_orr[1]:.6f}, {d_orr[2]:.6f}, {d_orr[3]:.6f}"
                oer_entry = f"{pH:.4f}, {voltage:.4f}, {d_oer[0]:.6f}, {d_oer[1]:.6f}, {d_oer[2]:.6f}, {d_oer[3]:.6f}"
                orr_data.append(orr_entry)
                oer_data.append(oer_entry)
            except Exception as e:
                print(f"Error calculating at pH={pH}, U={voltage}: {e}")

    if orr_data:
        orr_file = os.path.join(save_dir, f"{prefix}_ORR.txt")
        with open(orr_file, 'w') as f:
            f.write(f"# MICROKINETICS Parameters for ORR ({voltage_scale} voltage)\n")
            f.write(f"# Format: pH, Potential vs {voltage_scale} (V), dG1, dG2, dG3, dG4\n")
            f.write("\n".join(orr_data))
        print(f"ORR parameters saved to: {orr_file}")

    if oer_data:
        oer_file = os.path.join(save_dir, f"{prefix}_OER.txt")
        with open(oer_file, 'w') as f:
            f.write(f"# MICROKINETICS Parameters for OER ({voltage_scale} voltage)\n")
            f.write(f"# Format: pH, Potential vs {voltage_scale} (V), dG1, dG2, dG3, dG4\n")
            f.write("\n".join(oer_data))
        print(f"OER parameters saved to: {oer_file}")


def plot_reaction_pathway(output_dir, step_conditions, C_list, U_pzc_list, E0_list, E_zpe_st,
                          reaction_type="ORR", include_O2=True, voltage_scale="RHE"):
    step_dir = create_output_directory(os.path.join(output_dir, f"reaction_pathway_{reaction_type}"))

    plt.figure(figsize=(10, 8))
    colors = plt.cm.magma(np.linspace(0.2, 0.8, len(step_conditions)))

    step_names = ["Step 1", "Step 2", "Step 3", "Step 4"]
    step_desc = {
        "ORR": [
            "* + O2 + H+ + e- → *OOH",
            "*OOH + H+ + e- → *O + H2O",
            "*O + H+ + e- → *OH",
            "*OH + H+ + e- → * + H2O"
        ],
        "OER": [
            "* + H2O → *OH + H+ + e-",
            "*OH → *O + H+ + e-",
            "*O + H2O → *OOH + H+ + e-",
            "*OOH → * + O2 + H+ + e-"
        ]
    }[reaction_type]

    labels = [
        ["* + O₂", "*OOH", "*O", "*OH", "* + H₂O"],
        ["* + H2O", "*OH", "*O", "*OOH", "* + O2"]
    ][0 if reaction_type == "ORR" else 1]

    all_rds_results = {}
    error_occurred = False

    for cond_idx, cond in enumerate(step_conditions):
        label = cond.get('label', f'Condition {cond_idx + 1}')
        pH = cond['pH']
        voltage = cond['U']

        try:
            print(f"\nProcessing pathway for: {label} (U={voltage} {voltage_scale}, pH={pH})")

            if voltage_scale == "RHE":
                voltage_she = U_rhe_to_U_she(voltage, pH)
                print(f"  Converted voltage: RHE {voltage:.4f} V → SHE {voltage_she:.4f} V (pH={pH})")
            else:
                voltage_she = voltage

            delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
                voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st,
                include_O2
            )
            print(
                f"  Calculated intermediates: ΔG*OOH={delta_GOOH:.4f} eV, ΔG*O={delta_GO:.4f} eV, ΔG*OH={delta_GOH:.4f} eV")

            dG_steps = calculate_steps_from_intermediates(
                delta_GOOH, delta_GO, delta_GOH,
                voltage,
                reaction_type
            )

            if reaction_type == "ORR":
                initial_height = delta_G_water - 4 * voltage
            else:
                initial_height = 0
            print(f"  Initial height: {initial_height:.4f} eV",
                  f"(Total reaction ΔG = {delta_G_water} eV, U = {voltage} V)")

            print(f"  Step energies: ΔG1={dG_steps[0]:.4f} eV, ΔG2={dG_steps[1]:.4f} eV, "
                  f"ΔG3={dG_steps[2]:.4f} eV, ΔG4={dG_steps[3]:.4f} eV")
            print(f"  Initial height for {reaction_type}: {initial_height:.4f} eV")

            cumulative_energy = [initial_height]
            current_energy = initial_height
            for dG in dG_steps:
                current_energy += dG
                cumulative_energy.append(current_energy)

            barriers = []
            for i in range(4):
                barrier = cumulative_energy[i + 1] - cumulative_energy[i]
                barriers.append(barrier)

            rds_idx = np.argmax(barriers)
            rds_step_name = step_names[rds_idx]
            rds_step_desc = step_desc[rds_idx]
            rds_barrier = barriers[rds_idx]
            rds_end_energy = cumulative_energy[rds_idx + 1]

            all_rds_results[label] = {
                "rds_step": rds_step_name,
                "rds_step_desc": rds_step_desc,
                "rds_idx": rds_idx,
                "rds_barrier": rds_barrier,
                "rds_end_energy": rds_end_energy,
                "voltage": float(voltage),
                "pH": float(pH),
                "cumulative_energy": [float(e) for e in cumulative_energy],
                "dG_steps": [float(d) for d in dG_steps]
            }

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
                color=colors[cond_idx],
                linewidth=2.5,
                label=f"{label} (RDS: {rds_step_name}, Barrier: {rds_barrier:.2f} eV)"
            )

            for i, energy in enumerate(cumulative_energy):
                plt.text(i, energy + 0.05, f"{energy:.2f}",
                         ha='center', fontsize=9, color=colors[cond_idx])

            print(f"  Pathway plotted for {label}")

        except Exception as e:
            print(f"ERROR: Calculation for {label} failed: {str(e)}")
            error_occurred = True
            continue

    if not error_occurred:
        plt.xticks(range(5), labels, fontsize=12)
        plt.xlim(-0.5, 4.5)
        plt.ylabel('Free Energy (eV)', fontsize=14)
        plt.title(f'{reaction_type} Reaction Pathway', fontsize=16)
        plt.grid(True, linestyle='--', alpha=0.7)

        plt.legend(fontsize=10, loc='best', framealpha=0.7)
        plt.tight_layout()

        fig_path = os.path.join(step_dir, f"{reaction_type}_pathway.pdf")
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.savefig(fig_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
        print(f"  ➤ Pathway diagram saved: {fig_path}")
        plt.close()

    rds_path = os.path.join(step_dir, f"{reaction_type}_rds_info.csv")
    with open(rds_path, 'w', encoding='utf-8') as f:
        f.write("Condition,U,pH,RDS_Step,RDS_Description,RDS_Barrier(eV),RDS_End_Energy(eV)\n")
        for label, data in all_rds_results.items():
            f.write(
                f"{label},{data['voltage']},{data['pH']},{data['rds_step']},\"{data['rds_step_desc']}\","
                f"{data['rds_barrier']},{data['rds_end_energy']}\n"
            )

    print(f"  ➤ RDS details saved: {rds_path}")

    full_data_path = os.path.join(step_dir, f"{reaction_type}_full_energy_data.csv")
    with open(full_data_path, 'w', encoding='utf-8') as f:
        f.write("Condition,Intermediate,FreeEnergy(eV)\n")
        for label, data in all_rds_results.items():
            intermediates = ["Initial", "*OOH", "*O", "*OH", "Final"]
            for i, energy in enumerate(data["cumulative_energy"]):
                f.write(f"{label},{intermediates[i]},{energy:.6f}\n")

    print(f"  ➤ Full energy data saved: {full_data_path}")

    if error_occurred:
        print("\n⚠ Warning: Some conditions failed to calculate, only available data has been saved")
    else:
        print(f"\n✓ {reaction_type} pathway analysis completed. Results saved to: {step_dir}")

    return all_rds_results

def calculate_onset_potential_simple_old1(pH, reaction_type, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2=True, voltage_scale="RHE"):
    step_desc = {
        "ORR": ["1: * + O₂ + H⁺ + e⁻ → *OOH", "2: *OOH + H⁺ + e⁻ → *O + H₂O",
                "3: *O + H⁺ + e⁻ → *OH", "4: *OH + H⁺ + e⁻ → * + H₂O"],
        "OER": ["1: * + H₂O → *OH + H⁺ + e⁻", "2: *OH → *O + H⁺ + e⁻",
                "3: *O + H₂O → *OOH + H⁺ + e⁻", "4: *OOH → * + O₂ + H⁺ + e⁻"]
    }

    ref_voltage = 0.0
    voltage_she = U_rhe_to_U_she(ref_voltage, pH)

    print(f"Determining RDS step for {reaction_type} at pH={pH}, U={ref_voltage} V (RHE)")

    try:
        delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
            voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st,
            include_O2
        )

        dG1, dG2, dG3, dG4 = calculate_steps_from_intermediates(
            delta_GOOH, delta_GO, delta_GOH, ref_voltage, reaction_type
        )

        step_values = [dG1, dG2, dG3, dG4]
        rds_value = max(step_values)
        rds_step = step_values.index(rds_value)
        rds_desc = step_desc[reaction_type][rds_step]

        print(f"  Step free energies (U=0 RHE):")
        print(f"    Step1: {dG1:.4f} eV, Step2: {dG2:.4f} eV")
        print(f"    Step3: {dG3:.4f} eV, Step4: {dG4:.4f} eV")
        print(f"  RDS determined: Step {rds_step + 1} ({rds_desc}), ΔG = {rds_value:.4f} eV")

    except Exception as e:
        print(f"Error calculating RDS step: {str(e)}")
        return None, None

    print(f"High-precision calculation: Finding {reaction_type} onset potential (making step {rds_step + 1} ΔG≈0)...")

    FINE_THRESHOLD = 0.005
    COARSE_STEPS = 1200
    FINE_STEPS = 1001

    if reaction_type == "ORR":
        coarse_range = np.linspace(3.0, 0.0, COARSE_STEPS)
    else:
        coarse_range = np.linspace(1.0, 3.0, COARSE_STEPS)

    fine_search_range = None
    prev_dG = None
    prev_voltage = None
    found_onset = None
    coarse_result = None
    coarse_dgs = []

    print("  Stage 1: Coarse scanning (0.0025V resolution)")
    for voltage_rhe in coarse_range:
        try:
            voltage_she = U_rhe_to_U_she(voltage_rhe, pH)
            delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
                voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2
            )
            step_values = calculate_steps_from_intermediates(
                delta_GOOH, delta_GO, delta_GOH, voltage_rhe, reaction_type
            )
            current_dG = step_values[rds_step]
            coarse_dgs.append((voltage_rhe, current_dG))

            if abs(current_dG) < FINE_THRESHOLD:
                coarse_result = voltage_rhe
                print(f"    Direct solution found: U={voltage_rhe:.5f} V | ΔG={current_dG:.6f} eV")
                break

            if prev_dG is not None:
                if (reaction_type == "ORR" and prev_dG >= 0 and current_dG <= 0) or \
                   (reaction_type == "OER" and prev_dG >= 0 and current_dG <= 0):
                    slope = (current_dG - prev_dG) / (voltage_rhe - prev_voltage)
                    # Removed slope filter for robustness
                    interpolated_U = voltage_rhe - current_dG / slope
                    if min(prev_voltage, voltage_rhe) <= interpolated_U <= max(prev_voltage, voltage_rhe):
                        coarse_result = interpolated_U
                        print(f"    Sign change detected between {prev_voltage:.4f}V and {voltage_rhe:.4f}V")
                        print(f"    Interpolation: U={interpolated_U:.5f} V")
                        break

            prev_dG = current_dG
            prev_voltage = voltage_rhe

        except Exception as e:
            print(f"    Calculation error at U={voltage_rhe:.4f}V: {str(e)}")
            continue

    if coarse_result is None:
        print("  No sign change detected; falling back to minimal |ΔG| in coarse grid")
        if coarse_dgs:
            min_idx = min(range(len(coarse_dgs)), key=lambda i: abs(coarse_dgs[i][1]))
            coarse_result = coarse_dgs[min_idx][0]
            print(f"    Fallback: U={coarse_result:.5f} V | min ΔG={coarse_dgs[min_idx][1]:.6f} eV")

    if coarse_result is not None:
        print("  Coarse search completed! Preparing for fine scan...")
        fine_range = 0.05
        fine_start = max(0.0, coarse_result - fine_range)
        fine_end = min(3.0, coarse_result + fine_range)
        fine_search_range = np.linspace(fine_start, fine_end, FINE_STEPS)
        print(f"  Fine search focused on U={fine_start:.4f}-{fine_end:.4f} V (0.0001V steps)")

    if fine_search_range is not None:
        print("  Stage 2: Precision scanning (0.0001V resolution)")
        closest_voltage = None
        closest_dG = float('inf')

        for voltage_rhe in fine_search_range:
            try:
                voltage_she = U_rhe_to_U_she(voltage_rhe, pH)
                delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
                    voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2
                )
                step_values = calculate_steps_from_intermediates(
                    delta_GOOH, delta_GO, delta_GOH, voltage_rhe, reaction_type
                )
                current_dG = step_values[rds_step]
                abs_dG = abs(current_dG)

                if abs_dG < closest_dG:
                    closest_dG = abs_dG
                    closest_voltage = voltage_rhe

                if abs(current_dG) < FINE_THRESHOLD:
                    found_onset = voltage_rhe
                    print(f"    Precision solution found: U={voltage_rhe:.5f} V | ΔG={current_dG:.6f} eV")
                    break

                idx = list(fine_search_range).index(voltage_rhe)
                if idx % 100 == 0:
                    print(f"    Scan point {idx}/{len(fine_search_range)}: U={voltage_rhe:.5f}V | ΔG_rds={current_dG:.6f}eV")

            except Exception as e:
                continue

        if found_onset is None and closest_voltage is not None:
            found_onset = closest_voltage
            print(f"    Best solution: U={closest_voltage:.5f}V | min ΔG={closest_dG:.6f}eV")

    if found_onset is None:
        print("  Failed to find solution, using coarse result")
        found_onset = coarse_result

    if found_onset is None:
        print("  ERROR: Unable to determine onset potential")
    else:
        print(f"✓ Final onset potential: U={found_onset:.5f}V (RHE)")

    return rds_step, found_onset

def calculate_onset_potential_simple(pH, reaction_type, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2=True, voltage_scale="RHE"):
    step_desc = {
        "ORR": ["1: * + O₂ + H⁺ + e⁻ → *OOH", "2: *OOH + H⁺ + e⁻ → *O + H₂O",
                "3: *O + H⁺ + e⁻ → *OH", "4: *OH + H⁺ + e⁻ → * + H₂O"],
        "OER": ["1: * + H₂O → *OH + H⁺ + e⁻", "2: *OH → *O + H⁺ + e⁻",
                "3: *O + H₂O → *OOH + H⁺ + e⁻", "4: *OOH → * + O₂ + H⁺ + e⁻"]
    }
    ref_voltage = 0.0
    voltage_she = U_rhe_to_U_she(ref_voltage, pH)
    print(f"Determining RDS step for {reaction_type} at pH={pH}, U={ref_voltage} V (RHE)")
    try:
        delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
            voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st,
            include_O2
        )
        dG1, dG2, dG3, dG4 = calculate_steps_from_intermediates(
            delta_GOOH, delta_GO, delta_GOH, ref_voltage, reaction_type
        )
        step_values = [dG1, dG2, dG3, dG4]
        rds_value = max(step_values)
        rds_step = step_values.index(rds_value)
        rds_desc = step_desc[reaction_type][rds_step]
        print(f" Step free energies (U=0 RHE):")
        print(f" Step1: {dG1:.4f} eV, Step2: {dG2:.4f} eV")
        print(f" Step3: {dG3:.4f} eV, Step4: {dG4:.4f} eV")
        print(f" RDS determined: Step {rds_step + 1} ({rds_desc}), ΔG = {rds_value:.4f} eV")
    except Exception as e:
        print(f"Error calculating RDS step: {str(e)}")
        return None, None

    # 封装计算ΔG步骤
    def compute_dgs(U_rhe):
        voltage_she = U_rhe_to_U_she(U_rhe, pH)
        delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
            voltage_she, pH, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2
        )
        return calculate_steps_from_intermediates(delta_GOOH, delta_GO, delta_GOH, U_rhe, reaction_type)

    # 改进: 扩展自洽迭代 - 试所有可能RDS作为起始假设，选择最佳候选
    print(f"Self-consistent calculation: Finding {reaction_type} onset potential via RDS iteration...")
    candidates = []  # 收集 (U, max_dG, rds_idx, all_dgs)
    tolerance = 0.01  # ΔG ≈0 容差
    max_iters_per_rds = 5  # 每个RDS最多迭代
    for init_rds in range(4):  # 试每个步作为初始RDS
        current_rds = init_rds
        iter_count = 0
        converged = False
        while not converged and iter_count < max_iters_per_rds:
            iter_count += 1
            print(f"  Sub-iter {iter_count} for init RDS=Step {current_rds + 1}")
            # 二分法求当前RDS ΔG=0的U
            def target_for_rds(U):
                dgs = compute_dgs(U)
                return dgs[current_rds]
            try:
                # 调整范围: ORR [0,3] (max U), OER [-0.5,3.5] (min U)
                u_low, u_high = (0.0, 3.0) if reaction_type == "ORR" else (-0.5, 3.5)
                f_low = target_for_rds(u_low)
                f_high = target_for_rds(u_high)
                if f_low * f_high > 0:  # 同号，无根
                    break
                U_candidate = brentq(target_for_rds, u_low, u_high, xtol=1e-6)  # 精确根 (scipy)
                dgs_candidate = compute_dgs(U_candidate)
                dg_rds = dgs_candidate[current_rds]
                max_dg = np.max(dgs_candidate)
                if abs(dg_rds) > tolerance:
                    break
                if max_dg > tolerance:  # 其他步有垒
                    break
                new_rds = np.argmax(dgs_candidate)
                if new_rds == current_rds:
                    converged = True
                    candidates.append((U_candidate, max_dg, current_rds, dgs_candidate))
                    print(f"    ✓ Sub-converged for RDS {current_rds + 1}: U={U_candidate:.5f}V, max ΔG={max_dg:.6f}")
                else:
                    current_rds = new_rds
            except Exception as e:
                print(f"    Error in sub-iter: {e}")
                break

    # 选择最佳候选: 优先 max_dG <=0 的最小|U| (OER min U, ORR max U), 否则 min max_dG
    if candidates:
        valid_cands = [c for c in candidates if c[1] <= 0]
        if valid_cands:
            if reaction_type == "OER":
                best = min(valid_cands, key=lambda x: x[0])  # min U
            else:
                best = max(valid_cands, key=lambda x: x[0])  # max U
        else:
            best = min(candidates, key=lambda x: x[1])  # min max_dG (even if >0)
        found_onset, max_dg_best, final_rds_step, final_dgs = best
        final_rds_desc = step_desc[reaction_type][final_rds_step]
        print(f" ✓ Best from iteration: RDS = Step {final_rds_step + 1} ({final_rds_desc}) at U={found_onset:.5f} V, max ΔG={max_dg_best:.6f} eV")
        print(f"   All ΔG at onset: Step1={final_dgs[0]:.4f}, Step2={final_dgs[1]:.4f}, Step3={final_dgs[2]:.4f}, Step4={final_dgs[3]:.4f}")
        # 如果max >0.1, 警告物理限制
        if max_dg_best > 0.1:
            print(f"   WARNING: Persistent barrier {max_dg_best:.4f} eV - possible scaling relation limit or non-linear effects")
    else:
        print(" No candidates from iteration; fallback to optimized scan")
        found_onset = None
        final_rds_step = None

    # 改进Fallback: 用优化代替网格 - 最小化 max(ΔG(U))
    if found_onset is None:
        print(" Falling back to optimized search...")
        def obj(U):
            dgs = compute_dgs(U)
            penalty = 0.01 * (U - 1.23)**2 if reaction_type == "OER" else 0  # Penalize high overpotential for OER
            return np.max(dgs) + penalty

        u_bounds = (0.0, 3.0) if reaction_type == "ORR" else (-0.5, 4.0)  # 扩展 OER 范围
        try:
            f_low = obj(u_bounds[0])
            f_high = obj(u_bounds[1])
            if f_low * f_high < 0 or (f_low > 0 and f_high <= 0):
                def f_for_root(U):
                    return obj(U)
                found_onset = brentq(f_for_root, u_bounds[0], u_bounds[1], xtol=1e-6)
                print(f" Root found: U={found_onset:.5f}V where max ΔG≈0")
            else:
                res = minimize_scalar(obj, bounds=u_bounds, method='bounded', tol=1e-6)
                found_onset = res.x
                max_dg_best = np.max(compute_dgs(found_onset))  # Raw max, ignore penalty
                print(f" No zero-crossing; optimized U={found_onset:.5f}V, min max ΔG={max_dg_best:.6f}eV")
                if max_dg_best > 0.3:
                    print(f"   Scaling limit ~{max_dg_best:.3f} eV; fallback to linear U={rds_value:.3f} V")
                    found_onset = rds_value  # Linear approx: U = ΔG_RDS(0)
                    max_dg_best = 0.0
        except Exception as e:
            print(f" Optimization error: {e}; using initial RDS U")
            found_onset = rds_value  # 简单fallback: U = ΔG_RDS (线性假设)
            max_dg_best = 0.0  # 近似

        final_dgs = compute_dgs(found_onset)
        final_rds_step = np.argmax(final_dgs)
        final_rds_desc = step_desc[reaction_type][final_rds_step]
        print(f" Final RDS at onset U: Step {final_rds_step + 1} ({final_rds_desc})")
        print(f" All ΔG at onset: Step1={final_dgs[0]:.4f}, Step2={final_dgs[1]:.4f}, Step3={final_dgs[2]:.4f}, Step4={final_dgs[3]:.4f}")

    if found_onset is None:
        print(" ERROR: Unable to determine onset potential")
        return None, None
    else:
        print(f"✓ Final onset potential: U={found_onset:.5f}V (RHE)")
    return final_rds_step, found_onset

def U_rds_value(U, rds_step, pH, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2, voltage_scale, reaction_type):
    try:
        delta_GOOH, delta_GO, delta_GOH = calculate_intermediates(
            U, pH, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2
        )

        step_values = calculate_steps_from_intermediates(
            delta_GOOH, delta_GO, delta_GOH, U, reaction_type
        )

        return step_values[rds_step]
    except Exception as e:
        return float('inf') if reaction_type == "OER" else float('-inf')


def calculate_and_output_onset_potentials(output_dir, analysis_pH_values, C_list, U_pzc_list, E0_list, E_zpe_st, include_O2=True, voltage_scale="RHE"):
    print("\n" + "=" * 70)
    print("Calculating Onset Potentials")
    print("=" * 70)

    onset_dir = create_output_directory(os.path.join(output_dir, "onset_potentials"))

    onset_results = {
        "ORR": {},
        "OER": {},
        "metadata": {
            "voltage_scale": voltage_scale,
            "include_O2": include_O2
        }
    }

    csv_path = os.path.join(onset_dir, "onset_potentials.csv")

    with open(csv_path, 'w', encoding='utf-8') as f:
        f.write("Reaction,pH,RDS Step,RDS Description,Onset Potential (V vs RHE)\n")

    step_desc = {
        "ORR": ["1: * + O2 + H+ + e- → *OOH", "2: *OOH + H+ + e- → *O + H2O",
                "3: *O + H+ + e- → *OH", "4: *OH + H+ + e- → * + H2O"],
        "OER": ["1: * + H2O → *OH + H+ + e-", "2: *OH → *O + H+ + e-",
                "3: *O + H2O → *OOH + H+ + e-", "4: *OOH → * + O2 + H+ + e-"]
    }

    for pH in analysis_pH_values:
        print(f"\n> Calculating onset potentials for pH = {pH}")

        print(f"  ORR analysis:")
        orr_rds_step, orr_onset_potential = calculate_onset_potential_simple(
            pH, "ORR", C_list, U_pzc_list, E0_list, E_zpe_st,
            include_O2, voltage_scale
        )

        if orr_rds_step is not None:
            orr_rds_step = int(orr_rds_step)
            orr_desc = step_desc["ORR"][orr_rds_step]
        else:
            orr_desc = "Unknown"

        print(f"  OER analysis:")
        oer_rds_step, oer_onset_potential = calculate_onset_potential_simple(
            pH, "OER", C_list, U_pzc_list, E0_list, E_zpe_st,
            include_O2, voltage_scale
        )

        if oer_rds_step is not None:
            oer_rds_step = int(oer_rds_step)
            oer_desc = step_desc["OER"][oer_rds_step]
        else:
            oer_desc = "Unknown"

        onset_results["ORR"][str(pH)] = {
            "rds_step": orr_rds_step,
            "rds_desc": orr_desc if orr_rds_step is not None else None,
            "onset_potential": float(orr_onset_potential) if orr_onset_potential is not None else None
        }

        onset_results["OER"][str(pH)] = {
            "rds_step": oer_rds_step,
            "rds_desc": oer_desc if oer_rds_step is not None else None,
            "onset_potential": float(oer_onset_potential) if oer_onset_potential is not None else None
        }

        print(f"\nResults for pH={pH}:")
        if orr_onset_potential is not None:
            print(f"  ORR: RDS step = {orr_desc}, Onset potential = {orr_onset_potential:.4f} V")
        else:
            print(f"  ORR: Unable to determine onset potential")

        if oer_onset_potential is not None:
            print(f"  OER: RDS step = {oer_desc}, Onset potential = {oer_onset_potential:.4f} V")
        else:
            print(f"  OER: Unable to determine onset potential")

        orr_onset_str = f"{orr_onset_potential:.4f}" if orr_onset_potential is not None else "N/A"
        oer_onset_str = f"{oer_onset_potential:.4f}" if oer_onset_potential is not None else "N/A"

        with open(csv_path, 'a', encoding='utf-8') as f:
            f.write(f"ORR,{pH},{orr_rds_step},{orr_desc},{orr_onset_str}\n")
            f.write(f"OER,{pH},{oer_rds_step},{oer_desc},{oer_onset_str}\n")

    json_path = os.path.join(onset_dir, "onset_potentials.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        def convert(o):
            if isinstance(o, np.generic):
                return o.item()
            return o

        json.dump(onset_results, f, indent=4, default=convert)

    print(f"\nSaved onset potentials to: {onset_dir}")

    print("\nSummary of Onset Potentials:")
    print("pH | ORR Onset (V) | ORR RDS | OER Onset (V) | OER RDS")
    print("---------------------------------------------------------")
    for pH in analysis_pH_values:
        str_pH = str(pH)
        orr_result = onset_results["ORR"].get(str_pH, {})
        oer_result = onset_results["OER"].get(str_pH, {})

        orr_onset = orr_result.get("onset_potential")
        orr_rds = orr_result.get("rds_desc", "?")
        if not orr_rds: orr_rds = "?"

        oer_onset = oer_result.get("onset_potential")
        oer_rds = oer_result.get("rds_desc", "?")
        if not oer_rds: oer_rds = "?"

        orr_display = f"{orr_onset:.4f}" if orr_onset is not None else "N/A"
        oer_display = f"{oer_onset:.4f}" if oer_onset is not None else "N/A"

        orr_rds_short = orr_rds.split(':')[1].strip()[:15] if ':' in orr_rds else orr_rds[:15]
        oer_rds_short = oer_rds.split(':')[1].strip()[:15] if ':' in oer_rds else oer_rds[:15]

        print(
            f"{pH:2d} | {orr_display:>12} | {orr_rds_short}{'...' if len(orr_rds_short) == 15 else '':<3} | {oer_display:>12} | {oer_rds_short}{'...' if len(oer_rds_short) == 15 else '':<3}")

    return onset_results

def save_single_pH_results(output_dir, results, pH, reaction_type, voltage_scale):
    if not results or results["voltage"].size == 0:  # Fixed: Check if results is None or voltage array is empty
        print(f"Warning: No valid results for pH={pH}, {reaction_type}. Skipping save.")
        return

    filename = f"{reaction_type}_analysis_pH{pH}_VS_{voltage_scale}.csv"
    path = os.path.join(output_dir, filename)
    data = {
        f"Voltage ({voltage_scale})": results["voltage"],
        "delta_GOOH": results["intermediates"]["delta_GOOH"],
        "delta_GO": results["intermediates"]["delta_GO"],
        "delta_GOH": results["intermediates"]["delta_GOH"],
        "d_G1": results["steps"]["d_G1"],
        "d_G2": results["steps"]["d_G2"],
        "d_G3": results["steps"]["d_G3"],
        "d_G4": results["steps"]["d_G4"]
    }
    df = pd.DataFrame(data)
    df.to_csv(path, index=False)
    print(f"Saved {reaction_type} analysis data: {path}")