import numpy as np
import pandas as pd


def run_selected_games_analysis(inhibition_matrices, ic50_map, beta=0.05):
    """Runs the game-theoretic analysis for the 6 hand-selected key drug combinations

    across cell lines (A375 and HT29), computes payoff matrices using asymmetric
    game theory,
    identifies Pure Nash Equilibria (or Pareto Optimal fallback when absent),
    and returns a comprehensive summary dataframe.
    """
    # Define the 6 hand-selected key games based on biological relevance and literature
    selected_games = [
        {
            "key": "Lapatinib__BEZ-235__HT29",
            "label": "2x3",
            "transpose": False,
        },
        {
            "key": "Etoposide__Dasatinib__HT29",
            "label": "3x4",
            "transpose": False,
        },
        {
            "key": "5-FU__BEZ-235__HT29",
            "label": "4x4",
            "transpose": False,
        },
        {
            "key": "Dasatinib__BEZ-235__A375",
            "label": "2x2",
            "transpose": False,
        },
        {
            "key": "AZD1775__geldanamycin__A375",
            "label": "3x3",
            "transpose": False,
        },
        {
            "key": "Temozolomide__BEZ-235__A375",
            "label": "2x4",
            "transpose": True,  # Transpose to maintain N x M format (smaller dimension first)
        },
    ]

    print("=" * 80)
    print("          ASYMMETRIC GAME THEORY & NASH EQUILIBRIUM (SELECTED GAMES)")
    print("=" * 80 + "\n")

    game_summary_results = []

    # Iterate through each selected game
    for item in selected_games:
        raw_key = item["key"]
        label = item["label"]
        is_transposed = item["transpose"]

        # Retrieve raw inhibition matrix and split keys
        raw_matrix = inhibition_matrices[raw_key]
        dA_orig, dB_orig, cell = raw_key.split("__")

        # Handle matrix transposition if necessary to keep dimensions consistent
        if is_transposed:
            matrix_inhibition = raw_matrix.T
            dA, dB = dB_orig, dA_orig
        else:
            matrix_inhibition = raw_matrix
            dA, dB = dA_orig, dB_orig

        # Fetch reference IC50 values for both drugs on the specific cell line
        ic50_A = ic50_map[(dA, cell)]
        ic50_B = ic50_map[(dB, cell)]

        rows_doses_A = matrix_inhibition.index.values
        cols_doses_B = matrix_inhibition.columns.values
        rows, cols = matrix_inhibition.shape

        # Initialize payoff matrices
        payoff_A = np.zeros((rows, cols))
        payoff_B = np.zeros((rows, cols))

        # Compute Payoff values considering drug inhibition and penalty factor (beta)
        for r, dA_val in enumerate(rows_doses_A):
            for c, dB_val in enumerate(cols_doses_B):
                inh = matrix_inhibition.iloc[r, c]
                payoff_A[r, c] = inh - beta * (dA_val / ic50_A)
                payoff_B[r, c] = inh - beta * (dB_val / ic50_B)

        # Extract Pure Nash Equilibrium candidates
        nash_candidates = []
        for r in range(rows):
            for c in range(cols):
                is_A_best = payoff_A[r, c] == np.max(payoff_A[:, c])
                is_B_best = payoff_B[r, c] == np.max(payoff_B[r, :])

                if is_A_best and is_B_best:
                    nash_candidates.append(
                        {
                            "row_idx": r,
                            "col_idx": c,
                            "dose_A": rows_doses_A[r],
                            "dose_B": cols_doses_B[c],
                            "payoff_A": payoff_A[r, c],
                            "payoff_B": payoff_B[r, c],
                            "total_payoff": payoff_A[r, c] + payoff_B[r, c],
                        }
                    )

        # Handle equilibrium selection (Pure Nash vs Pareto Optimal fallback if no pure Nash exists)
        if len(nash_candidates) > 0:
            sorted_candidates = sorted(
                nash_candidates,
                key=lambda x: (-x["total_payoff"], x["dose_A"] + x["dose_B"]),
            )
            best_nash = sorted_candidates[0]
            eq_type = (
                "Unique Nash" if len(nash_candidates) == 1 else "Pareto Optimal Nash"
            )
        else:
            # Fallback to max total payoff (Pareto Optimal cell) when no pure Nash equilibrium is found
            all_cells = []
            for r in range(rows):
                for c in range(cols):
                    all_cells.append(
                        {
                            "row_idx": r,
                            "col_idx": c,
                            "dose_A": rows_doses_A[r],
                            "dose_B": cols_doses_B[c],
                            "payoff_A": payoff_A[r, c],
                            "payoff_B": payoff_B[r, c],
                            "total_payoff": payoff_A[r, c] + payoff_B[r, c],
                        }
                    )
            sorted_cells = sorted(
                all_cells,
                key=lambda x: (-x["total_payoff"], x["dose_A"] + x["dose_B"]),
            )
            best_nash = sorted_cells[0]
            eq_type = "Pareto Optimal (No Pure Nash)"

        # Construct Payoff Matrix DataFrame for terminal printing
        payoff_matrix_df = pd.DataFrame(
            index=[f"Dose {dA}: {d}" for d in rows_doses_A],
            columns=[f"Dose {dB}: {d}" for d in cols_doses_B],
        )

        for r in range(rows):
            for c in range(cols):
                val_a = round(payoff_A[r, c], 4)
                val_b = round(payoff_B[r, c], 4)
                is_nash = (r == best_nash["row_idx"]) and (c == best_nash["col_idx"])
                mark = " ★ " if is_nash else ""
                payoff_matrix_df.iloc[r, c] = f"({val_a}, {val_b}){mark}"

        # Print results for the current game
        print(f"🎮 [Game {label}] -> {dA} + {dB} | Cell Line: {cell}")
        print(
            f"   Matrix Dimensions: {rows}x{cols} | IC50 {dA}: {ic50_A} µM | IC50"
            f" {dB}: {ic50_B} µM"
        )
        print(
            f'   Optimal Strategy: Dose {dA} = {best_nash["dose_A"]} µM, Dose {dB}'
            f' = {best_nash["dose_B"]} µM'
        )
        print(
            f'   Payoff {dA} = {round(payoff_nash_A := best_nash["payoff_A"], 4)},'
            f' Payoff {dB} = {round(best_nash["payoff_B"], 4)} | Total ='
            f' {round(best_nash["total_payoff"], 4)}'
        )
        print("\n Payoff Matrix Table (U_A, U_B):")
        print(payoff_matrix_df.to_string())
        print("-" * 80 + "\n")

        # Append summary metrics for final table
        game_summary_results.append(
            {
                "Size": label,
                "Drug A": dA,
                "Drug B": dB,
                "Cell Line": cell,
                "Opt Dose A (µM)": best_nash["dose_A"],
                "Opt Dose B (µM)": best_nash["dose_B"],
                "Payoff A": round(best_nash["payoff_A"], 4),
                "Payoff B": round(best_nash["payoff_B"], 4),
                "Total Payoff": round(best_nash["total_payoff"], 4),
                "Equilibrium": eq_type,
            }
        )

    summary_df = pd.DataFrame(game_summary_results)
    print(" Summary Table of 6 Selected Games for Manuscript:")
    print(summary_df.to_string())
    return summary_df
