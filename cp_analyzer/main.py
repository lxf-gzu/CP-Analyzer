# main.py
import argparse
import sys
from pathlib import Path
import os
import logging
import time
import textwrap
import importlib.util

# ====================== CP-Analyzer 初始化 ======================
from cp_analyzer import initialize_config

# ====================== 包内导入 ======================
from cp_analyzer.tasks import TASK_REGISTRY
from cp_analyzer.utils import *
from cp_analyzer.data_processing import *
from cp_analyzer.calculations import *
from cp_analyzer.plotting import *


def load_user_config():
    """动态加载当前目录下的 config.py"""
    cwd = Path.cwd()
    config_path = cwd / "config.py"

    if not config_path.exists():
        initialize_config()

    if not config_path.exists():
        print("❌ Failed to create config.py")
        sys.exit(1)

    spec = importlib.util.spec_from_file_location("user_config", config_path)
    user_config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(user_config)

    globals().update({k: v for k, v in user_config.__dict__.items() if not k.startswith("__")})
    print(f"✅ Loaded config from: {config_path}")
    return user_config

def main():
    load_user_config()

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

    def check_task_dependencies(selected_tasks, task_dependencies):
        selected_tasks = set(selected_tasks)

        for task_id in sorted(selected_tasks):
            if task_id in task_dependencies:
                required = set(task_dependencies[task_id])
                missing = required - selected_tasks

                if missing:
                    missing_str = " ".join(str(t) for t in sorted(missing))
                    required_str = " ".join(str(t) for t in sorted(required))
                    selected_str = " ".join(str(t) for t in sorted(selected_tasks))

                    print("\n" + "=" * 80)
                    print(f"ERROR: TASK {task_id} depends on TASK {required_str}")
                    print(f"You selected: --task {selected_str}")
                    print("Please run:")
                    print(f"  python main.py --task {missing_str} {task_id}")
                    print("=" * 80 + "\n")

                    logging.error(
                        f"Task dependency error: Task {task_id} requires {required_str}, "
                        f"but selected tasks are {selected_str}"
                    )
                    return False

        return True

    # Detailed task descriptions for help message
    task_descriptions = {
        1: """Task 1: Obtaining Fitting Parameters
           Purpose: Fits electrochemical data to obtain coefficients and parameters.
           Usage: Specify --data-source 'file' or 'database' and --file-paths for input files.
           Config: Set DEFAULT_FILE_PATHS, DEFAULT_DATA_SOURCE in config.py.
           Output: Fitting parameters saved in output_dir.""",
        2: """Task 2: Plotting Free Energy vs SHE Potential
           Purpose: Plots free energy as a function of SHE potential for electrochemical systems.
           Usage: Requires Task 1 data and --data-source 'file'.
           Config: Set DEFAULT_FILE_PATHS, DEFAULT_VOLTAGE_SCALE in config.py.
           Output: Free energy plots and data saved in output_dir/free_energy_fits.""",
        3: """Task 3: Q vs. pH Analysis
           Purpose: Analyzes charge (Q) as a function of pH for electrochemical systems.
           Usage: Enable with --run-q-vs-ph. Requires --data-source 'file' and valid --file-paths.
           Config: Set DEFAULT_FILE_PATHS in config.py for input files.
           Output: Q vs. pH plots/data saved in output_dir/Q_vs_pH.""",
        4: """Task 4: Onset Potential Analysis
           Purpose: Calculates onset potentials for electrochemical reactions across specified pH values.
           Usage: Enable with --run-onset-analysis. Requires Task 1 data and --voltage-scale 'RHE'.
           Config: Set DEFAULT_ONSET_PH_VALUES, DEFAULT_E_ZPE_ST_LITERATURE/DEFAULT_E_ZPE_ST_FILE in config.py.
           Output: Onset potential data/plots saved in output_dir/onset_potentials.""",
        5: """Task 5: Single pH Point Analysis
           Purpose: Analyzes ORR and OER reaction energetics at specific pH values.
           Usage: Requires Task 1 data, --data-source 'file', and --run-onset-analysis.
           Config: Set DEFAULT_ONSET_PH_VALUES, DEFAULT_VOLTAGE_SCALE in config.py.
           Output: ORR/OER energy data (csv) saved in output_dir/pH_analysis.""",
        6: """Task 6: Heatmap Analysis
           Purpose: Generates heatmaps of ORR and OER properties across pH and potential.
           Usage: Requires Task 1 data and --data-source 'file'.
           Config: Set DEFAULT_PH_RANGE, DEFAULT_VOLTAGE_RANGE in config.py.
           Output: ORR/OER heatmap plots (png) saved in output_dir/heatmaps.""",
        7: """Task 7: Reaction Pathway Diagram Analysis
           Purpose: Generates ORR and OER reaction pathway diagrams at specific pH and voltage points.
           Usage: Requires Task 1 data and --data-source 'file'.
           Config: Set DEFAULT_ONSET_PH_VALUES, DEFAULT_VOLTAGE_SCALE in config.py.
           Output: Pathway diagrams (pdf) and energy data (csv) saved in output_dir/reaction_pathway_{ORR,OER}.""",
        8: """Task 8: Microkinetics Parameter Export
           Purpose: Exports microkinetic parameters for ORR and OER reactions.
           Usage: Requires Task 1 data and --data-source 'file'.
           Config: Set DEFAULT_ONSET_PH_VALUES in config.py.
           Output: Microkinetic parameters (txt) saved in output_dir/microkinetics.""",
        9: """Task 9: Charge vs Potential Analysis
           Purpose: Analyzes charge (Q) as a function of potential for electrochemical systems, providing linear fit equations and VASP NELECT adjustment guidance.
           Usage: Requires Task 1 data and --data-source 'file'. Specify voltage scale with --voltage-scale (default: RHE).
           Config: Set DEFAULT_VOLTAGE_SCALE, DEFAULT_PH, DEFAULT_FILE_PATHS in config.py.
           Output: Charge vs. potential data (txt), fitting results (csv), combined plot (png), and VASP NELECT instructions saved in output_dir/charge_analysis.""",
        10: """Task 10: Generic Adsorption Systems Analysis
            Purpose: Analyzes adsorption energies and fitting parameters for generic systems (e.g., *CO).
            Usage: Requires --data-source 'file' and valid --file-paths.
            Config: Set DEFAULT_FILE_PATHS in config.py.
            Output: Adsorption fitting parameters (txt/csv) saved in output_dir.""",
        11: """Task 11: Create VASP Calculation Directories
            Purpose: Creates directories for VASP calculations with varying NELECT values.
            Usage: Requires neutral system NELECT value.
            Config: Set DEFAULT_NELECT in config.py.
            Output: VASP directories created in output_dir/vasp_calculations.""",
        12: """Task 12: Extract VASP Data to Excel
            Purpose: Extracts data from VASP LOCPOT and OUTCAR files to Excel.
            Usage: Requires VASP output files in calculation directories.
            Config: Set DEFAULT_VASP_OUTPUT_DIR in config.py.
            Output: Excel file with extracted data saved in output_dir, or error if files are missing.""",
        13: """Task 13: ORR/OER Microkinetic Simulation
            Purpose: Performs microkinetic simulations for ORR and OER, calculating current density and coverages.
            Usage: Requires Task 1 data and --data-source 'file'.
            Config: Set DEFAULT_ONSET_PH_VALUES, DEFAULT_VOLTAGE_RANGE in config.py.
            Output: Current density, coverage data, and plots (png) saved in output_dir/Task13.""",
        14: """Task 14: OER/ORR Free Energy Diagram (Constant Charge)
            Purpose: Generates free energy diagrams for ORR and OER reactions at constant charge, computing intermediates, step energies with pH/U corrections, identifying RDS, and plotting reaction pathways.
            Usage: Enable with --run_oer_orr. Specify --pH, --U, and --reaction-type (ORR/OER). Requires VASP energies and corrections from config.
            Config: Set OER_ORR_CONFIG (pH, U, reaction_type, vasp_energies, corrections, ZPE_TDS, delta_G_water) in config.py; args override pH/U/reaction_type.
            Output: Pathway diagrams (pdf/png), intermediates/steps data (csv), RDS info (csv) saved in output_dir/oer_orr_analysis.""",
        15: """Task 15: OER/ORR Activity Volcano Plot
            Purpose: Generates activity volcano plots for OER and ORR.
            Usage: Requires Task 1 data and --data-source 'file'.
            Config: Set DEFAULT_PLOT_RANGE, DEFAULT_COLOR_RANGE in config.py.
            Output: Volcano plots (png) and data (csv) saved in output_dir/Task13/2D_Volcano."""
    }

    # Custom help formatter to include task descriptions
    class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
        def add_usage(self, usage, actions, groups, prefix=None):
            if prefix is None:
                prefix = "Usage: "
            super().add_usage(usage, actions, groups, prefix)

    # Initialize parser with custom description
    parser = argparse.ArgumentParser(
        description=textwrap.dedent("""\
            Electrochemical Free Energy Calculator
            -------------------------------------
            This script performs electrochemical free energy calculations and analyses for reactions like HER, OER, and ORR.
            Select tasks using --task (e.g., --task 1 2) or enable specific tasks with flags (e.g., --run-orr-oer-microkinetics).
            Note: --task N and --run-<task> flags (e.g., --task 13 and --run-orr-oer-microkinetics) are equivalent for tasks with dedicated flags; use either to enable. Dedicated flags provide additional control but are redundant with --task.
            Configure defaults in config.py (e.g., DEFAULT_DATA_SOURCE, DEFAULT_ANALYSIS_PH_VALUES).

            Available Tasks:
            ----------------
            {}
            """.format("\n".join([f"Task {k}: {v}" for k, v in task_descriptions.items()]))),
        formatter_class=CustomHelpFormatter
    )

    # Add arguments
    parser.add_argument('--output-dir', default=DEFAULT_OUTPUT_DIR, help="Output directory for results and logs")
    parser.add_argument('--data-source', default=DEFAULT_DATA_SOURCE, choices=["literature", "file", "literature_params"],
                        help="Data source: 'literature' (predefined), 'file' (Excel/CSV), or 'literature_params' (literature with custom params)")
    parser.add_argument('--voltage-scale', default=DEFAULT_VOLTAGE_SCALE, choices=["SHE", "RHE"],
                        help="Voltage reference: SHE (Standard Hydrogen Electrode) or RHE (Reversible Hydrogen Electrode)")
    parser.add_argument('--include-o2', action='store_true', default=DEFAULT_INCLUDE_O2,
                        help="Include O2-related calculations (e.g., for ORR)")
    parser.add_argument('--file-paths', nargs='+', default=DEFAULT_FILE_PATHS,
                        help="Excel/CSV files for 'file' data source")
    parser.add_argument('--task', nargs='+', type=int, default=list(TASK_REGISTRY.keys()),
                        help="Task IDs to run (e.g., 1 2). See task descriptions above.")
    parser.add_argument('--analysis-ph-values', nargs='+', type=float, default=DEFAULT_ANALYSIS_PH_VALUES,
                        help="pH values for analysis (e.g., 0 7 14)")
    parser.add_argument('--onset-ph-values', nargs='+', type=float, default=DEFAULT_ONSET_PH_VALUES,
                        help="pH values for onset potential analysis")
    parser.add_argument('--default-potential-points', nargs='+', type=float, default=DEFAULT_DEFAULT_POTENTIAL_POINTS,
                        help="Potential points for analysis (V vs. SHE/RHE)")
    parser.add_argument('--charge-ph', type=int, default=DEFAULT_CHARGE_PH,
                        help="pH value for charge-related calculations")
    parser.add_argument('--u-min', type=int, default=DEFAULT_U_MIN,
                        help="Minimum potential for analysis (V)")
    parser.add_argument('--u-max', type=int, default=DEFAULT_U_MAX,
                        help="Maximum potential for analysis (V)")
    parser.add_argument('--heatmap-ph-range', type=float, nargs=2, default=[0, 14],
                        help="Min and max pH for heatmaps (e.g., 0 14)")
    parser.add_argument('--heatmap-voltage-range', type=float, nargs=2, default=[0, 1.5],
                        help="Min and max voltage for heatmaps (e.g., 0 1.5)")
    parser.add_argument('--run-q-vs-ph', action='store_true', default=DEFAULT_RUN_Q_VS_PH,
                        help="Enable Task 4: Q vs. pH analysis")
    parser.add_argument('--run-single-ph', action='store_true', default=DEFAULT_RUN_SINGLE_PH,
                        help="Enable Task 5: Single pH free energy analysis")
    parser.add_argument('--run-heatmaps', action='store_true', default=DEFAULT_RUN_HEATMAPS,
                        help="Enable Task 6: pH-voltage heatmap generation")
    parser.add_argument('--run-microkinetics', action='store_true', default=DEFAULT_RUN_MICROKINETICS,
                        help="Enable Task 7: Microkinetic modeling")
    parser.add_argument('--run-step-diagrams', action='store_true', default=DEFAULT_RUN_STEP_DIAGRAMS,
                        help="Enable Task 8: Reaction step diagrams")
    parser.add_argument('--run-onset-analysis', action='store_true', default=DEFAULT_RUN_ONSET_ANALYSIS,
                        help="Enable Task 9: Onset potential analysis")
    parser.add_argument('--run-qvspotential', action='store_true', default=DEFAULT_RUN_QVSPOTENTIAL,
                        help="Enable Task 10: Q vs. potential analysis")
    parser.add_argument('--run-other-ads-data', action='store_true', default=DEFAULT_RUN_OTHER_ADS_DATA,
                        help="Enable processing of additional adsorption data")
    parser.add_argument('--run-vasp-dir-creation', action='store_true', default=DEFAULT_RUN_VASP_DIR_CREATION,
                        help="Enable Task 11: Create VASP directories with modified NELECT")
    parser.add_argument('--vasp-input-files', nargs='+', default=DEFAULT_VASP_INPUT_FILES,
                        help="VASP input files (default: INCAR KPOINTS POSCAR POTCAR jy-sub_54_.vasp)")
    parser.add_argument('--vasp-submit-script', default=DEFAULT_VASP_SUBMIT_SCRIPT,
                        help="VASP job submission script (default: jy-sub_54_.vasp)")
    parser.add_argument('--submit-vasp-jobs', action='store_true', default=DEFAULT_SUBMIT_VASP_JOBS,
                        help="Submit VASP jobs via sbatch after directory creation")
    parser.add_argument('--run-vasp-data-extraction', action='store_true', default=DEFAULT_RUN_VASP_DATA_EXTRACTION,
                        help="Enable Task 12: Extract VASP data to Excel")
    parser.add_argument('--run-orr-oer-microkinetics', action='store_true', default=DEFAULT_RUN_ORR_OER_MICROKINETICS,
                        help="Enable Task 13: ORR/OER microkinetic simulation")
    parser.add_argument('--run-oer-orr', action='store_true', default=False,
                        help="Enable Task 14: OER/ORR free energy diagram (requires --pH, --reaction-type, --U)")
    parser.add_argument('--pH', type=float, default=None,
                        help="pH value for Task 14 (OER/ORR diagram)")
    parser.add_argument('--reaction-type', type=str, choices=["ORR", "OER"], default=None,
                        help="Reaction type for Task 14 (ORR or OER)")
    parser.add_argument('--U', type=float, default=None,
                        help="Potential (V vs. SHE) for Task 14")
    parser.add_argument('--run-oer-orr-volcano', action='store_true', default=DEFAULT_RUN_OER_ORR_VOLCANO,
                        help="Enable Task 15: OER/ORR activity volcano plot")

    args = parser.parse_args()
    # ===============================
    # Task dependency validation
    # ===============================
    # Ensure output directory
    args.output_dir = os.path.abspath(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    if not check_task_dependencies(args.task, TASK_DEPENDENCIES):
        sys.exit(1)

    # Ensure output directory is created and normalized
    args.output_dir = os.path.abspath(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)
    args.heatmap_ph_range = np.linspace(args.heatmap_ph_range[0], args.heatmap_ph_range[1], 150)
    args.heatmap_voltage_range = np.linspace(args.heatmap_voltage_range[0], args.heatmap_voltage_range[1], 150)

    # Set up logging
    log_file = os.path.join(args.output_dir, 'execution_log.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info("Starting Electrochemical Free Energy Calculator")

    # Task status tracking
    task_status = {task_id: "Not Started" for task_id in TASK_REGISTRY.keys()}



    # Add logging for Task 14 parameters
    if 14 in args.task:
        logging.info("-" * 80)
        logging.info("Task 14 Parameters:".center(80))
        logging.info(
            f"{'pH:':<20} {args.pH if args.pH is not None else 'Using default pH ' + str(OER_ORR_CONFIG['pH'])}")
        logging.info(
            f"{'Potential:':<20} {args.U if args.U is not None else 'Using default U ' + str(OER_ORR_CONFIG['U']) + ' V'}")
        logging.info(
            f"{'Reaction type:':<20} {args.reaction_type if args.reaction_type else 'Using default reaction ' + OER_ORR_CONFIG['reaction_type']}")
        logging.info("-" * 80)

    logging.info("=" * 80)
    logging.info("ELECTROCHEMICAL FREE ENERGY CALCULATOR".center(80))
    logging.info("=" * 80)
    logging.info(f"{'Data Source:':<20} {args.data_source}")
    logging.info(f"{'Voltage Reference:':<20} {args.voltage_scale}")
    logging.info(f"{'Include O2:':<20} {'Yes' if args.include_o2 else 'No'}")
    logging.info(f"{'Output Directory:':<20} {args.output_dir}")
    logging.info("-" * 80)

    C_list = U_pzc_list = E0_list = R2_list = data_to_use = None
    E_zpe_st = DEFAULT_E_ZPE_ST_LITERATURE if args.data_source in ["literature", "literature_params"] else DEFAULT_E_ZPE_ST_FILE

    for task_id in sorted(args.task):
        if task_id not in TASK_REGISTRY:
            logging.warning(f"Invalid task {task_id}, skipping")
            task_status[task_id] = "Skipped"
            continue

        func = TASK_REGISTRY[task_id]
        start_time = time.time()
        logging.info(f"Executing Task {task_id}")
        try:
            if task_id == 1:
                C_list, U_pzc_list, E0_list, R2_list, data_to_use = func(args)
                task_status[task_id] = "Completed"
            elif task_id == 2:
                func(args, C_list, U_pzc_list, E0_list, R2_list, data_to_use)
                task_status[task_id] = "Completed"
            elif task_id == 3:
                func(args, args.file_paths)
                task_status[task_id] = "Completed"
            elif task_id == 4:
                func(args, C_list, U_pzc_list, E0_list, E_zpe_st)
                task_status[task_id] = "Completed"
            elif task_id == 5:
                func(args, C_list, U_pzc_list, E0_list, E_zpe_st)
                task_status[task_id] = "Completed"
            elif task_id == 6:
                func(args, C_list, U_pzc_list, E0_list, E_zpe_st)
                task_status[task_id] = "Completed"
            elif task_id == 7:
                func(args, C_list, U_pzc_list, E0_list, E_zpe_st)
                task_status[task_id] = "Completed"
            elif task_id == 8:
                func(args, C_list, U_pzc_list, E0_list, E_zpe_st)
                task_status[task_id] = "Completed"
            elif task_id == 9:
                func(args, args.file_paths)
                task_status[task_id] = "Completed"
            elif task_id == 10:
                func(args, args.file_paths)
                task_status[task_id] = "Completed"
            elif task_id == 11 and args.run_vasp_dir_creation:
                TASK_REGISTRY[task_id](args, args.vasp_input_files)
                task_status[task_id] = "Completed"
            elif task_id == 12 and args.run_vasp_data_extraction:
                TASK_REGISTRY[task_id](args)
                task_status[task_id] = "Completed"
            elif task_id == 13 and args.run_orr_oer_microkinetics:
                func(args)
                task_status[task_id] = "Completed"
                # Check Task 13 output
                expected_files = ['U_coverage_current_data.txt', 'Key_potentials.txt',
                                'O2_dl_Coverage.png', 'Intermediate_Coverages.png',
                                'ORR_OER_Polarization_Curve.png', 'TOF_Plot.png']
                for f in expected_files:
                    file_path = os.path.join(args.output_dir, f)
                    if not os.path.exists(file_path):
                        logging.warning(f"Expected output {f} not found in {args.output_dir}")

            elif task_id == 14:
                logging.info(f"Executing Task 14 (OER/ORR Free Energy Diagram) - "
                             f"Reaction type: {args.reaction_type if args.reaction_type else OER_ORR_CONFIG['reaction_type']}, "
                             f"pH: {args.pH if args.pH is not None else OER_ORR_CONFIG['pH']}, "
                             f"U: {args.U if args.U is not None else OER_ORR_CONFIG['U']} V")
                func(args)
                task_status[task_id] = "Completed"
            elif task_id == 15 and args.run_oer_orr_volcano:
                func(args)
                task_status[task_id] = "Completed"
        except Exception as e:
            logging.error(f"Task {task_id} failed: {e}")
            task_status[task_id] = f"Failed: {e}"

        execution_time = time.time() - start_time
        task_status[task_id] += f" (Time: {execution_time:.2f}s)"
        logging.info(f"Task {task_id} completed in {execution_time:.2f} seconds")

    logging.info("\n" + "=" * 80)
    logging.info("ALL ANALYSES COMPLETED SUCCESSFULLY!".center(80))
    logging.info("=" * 80)
    logging.info(f"Results saved to: {os.path.abspath(args.output_dir)}")
    logging.info("=" * 80)

    # Console output summary
    print("\n=== Task Execution Summary ===")
    for task_id, status in task_status.items():
        print(f"Task {task_id}: {status}")

if __name__ == "__main__":
    main()