import numpy as np
import pandas as pd


def solve_common_payoff_games(inhibition_matrices):
  """Model common-payoff games and find Nash equilibria based on inhibition matrices."""
  common_payoff_results = []

  for key, matrix_inhibition in inhibition_matrices.items():
    dA, dB, cell = key.split('__')
    val_matrix = matrix_inhibition.values
    rows, cols = val_matrix.shape

    if rows < 2 or cols < 2:
      continue

    nash_candidates = []
    for r in range(rows):
      for c in range(cols):
        payoff = val_matrix[r, c]
        is_row_best = payoff == np.max(val_matrix[:, c])
        is_col_best = payoff == np.max(val_matrix[r, :])

        if is_row_best and is_col_best:
          nash_candidates.append({
              'dose_A': matrix_inhibition.index[r],
              'dose_B': matrix_inhibition.columns[c],
              'payoff': payoff,
              'total_dose': matrix_inhibition.index[r]
              + matrix_inhibition.columns[c],
          })

    if nash_candidates:
      sorted_candidates = sorted(
          nash_candidates, key=lambda x: (-x['payoff'], x['total_dose'])
      )
      best_nash = sorted_candidates[0]

      common_payoff_results.append({
          'Drug_A': dA,
          'Drug_B': dB,
          'Cell_Line': cell,
          'Matrix_Shape': f'{rows}x{cols}',
          'Optimal_Dose_A_uM': best_nash['dose_A'],
          'Optimal_Dose_B_uM': best_nash['dose_B'],
          'Max_Inhibition_Payoff': round(best_nash['payoff'], 4),
          'Equilibrium_Type': (
              'Unique Nash'
              if len(nash_candidates) == 1
              else 'Pareto Optimal Nash'
          ),
      })

  return pd.DataFrame(common_payoff_results)


def solve_conflicting_payoff_games(inhibition_matrices, ic50_map, beta=0.05):
  """Model conflicting-payoff games with beta penalty parameter and find Nash equilibria."""
  conflicting_results = []

  for key, matrix_inhibition in inhibition_matrices.items():
    dA, dB, cell = key.split('__')
    rows_doses_A = matrix_inhibition.index.values
    cols_doses_B = matrix_inhibition.columns.values
    rows, cols = matrix_inhibition.shape

    if rows < 2 or cols < 2:
      continue

    ic50_A = ic50_map.get((dA, cell))
    ic50_B = ic50_map.get((dB, cell))
    if ic50_A is None or ic50_B is None:
      continue

    payoff_A = np.zeros(matrix_inhibition.shape)
    payoff_B = np.zeros(matrix_inhibition.shape)

    for r, dA_val in enumerate(rows_doses_A):
      for c, dB_val in enumerate(cols_doses_B):
        inhibition = matrix_inhibition.iloc[r, c]
        payoff_A[r, c] = inhibition - beta * (dA_val / ic50_A)
        payoff_B[r, c] = inhibition - beta * (dB_val / ic50_B)

    nash_candidates = []
    for r in range(len(rows_doses_A)):
      for c in range(len(cols_doses_B)):
        is_A_best = payoff_A[r, c] == np.max(payoff_A[:, c])
        is_B_best = payoff_B[r, c] == np.max(payoff_B[r, :])

        if is_A_best and is_B_best:
          nash_candidates.append({
              'dose_A': rows_doses_A[r],
              'dose_B': cols_doses_B[c],
              'payoff_A': payoff_A[r, c],
              'payoff_B': payoff_B[r, c],
              'total_payoff': payoff_A[r, c] + payoff_B[r, c],
          })

    if nash_candidates:
      sorted_candidates = sorted(
          nash_candidates,
          key=lambda x: (-x['total_payoff'], x['dose_A'] + x['dose_B']),
      )
      best_nash = sorted_candidates[0]

      conflicting_results.append({
          'Drug_A': dA,
          'Drug_B': dB,
          'Cell_Line': cell,
          'Matrix_Shape': f'{rows}x{cols}',
          'IC50_A_uM': ic50_A,
          'IC50_B_uM': ic50_B,
          'Optimal_Dose_A_uM': best_nash['dose_A'],
          'Optimal_Dose_B_uM': best_nash['dose_B'],
          'Payoff_Drug_A': round(best_nash['payoff_A'], 4),
          'Payoff_Drug_B': round(best_nash['payoff_B'], 4),
          'Total_Payoff': round(best_nash['total_payoff'], 4),
          'Equilibrium_Type': (
              'Unique Nash'
              if len(nash_candidates) == 1
              else 'Pareto Optimal Nash'
          ),
      })

  # Return dataframe after iterating through all matrices
  return pd.DataFrame(conflicting_results)