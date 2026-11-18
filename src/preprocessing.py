import warnings
import numpy as np
import pandas as pd

# Ignore warnings and set display options
warnings.filterwarnings('ignore')
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)


def load_datasets(synergy_path, single_agent_path):
  """Load pre-filtered synergy and single-agent datasets from processed data folder."""
  final_dataset = pd.read_csv("data/processed/filtered_ic50_synergy_data.csv")
  df_single = pd.read_csv("data/processed/filtered_single_agent_data.csv")
  return final_dataset, df_single


def calculate_ic50_map(df_single):
  """Calculate approximate empirical IC50 from single-agent data."""

  def calculate_ic50_nearest(group):
    group = group.copy()
    group['diff_from_50'] = np.abs(group['X/X0'] - 0.5)
    nearest_row = group.loc[group['diff_from_50'].idxmin()]
    closest_viability = nearest_row['X/X0']
    if 0.4 <= closest_viability <= 0.6:
      return nearest_row['Drug_concentration (µM)']
    else:
      return None

  # Ensure concentration column is numeric
  col_L_name = df_single.columns[11] if len(df_single.columns) > 11 else 'Drug_concentration (µM)'
  if col_L_name in df_single.columns:
    df_single[col_L_name] = pd.to_numeric(df_single[col_L_name], errors='coerce')

  sort_cols = [c for c in ['drug_name', 'cell_line', 'Drug_concentration (µM)'] if c in df_single.columns]
  if sort_cols:
    df_single = df_single.sort_values(by=sort_cols).reset_index(drop=True)

  # Create an IC50 mapping dictionary based on (drug_name, cell_line)
  ic50_map = (
      df_single.groupby(['drug_name', 'cell_line'])
      .apply(calculate_ic50_nearest)
      .to_dict()
  )
  return ic50_map


def build_inhibition_matrices(final_dataset, ic50_map):
  """Build inhibition matrices for unique drug combinations and cell lines."""
  unique_combos = final_dataset[
      ['drugA_name', 'drugB_name', 'cell_line']
  ].drop_duplicates()
  inhibition_matrices = {}

  for _, row in unique_combos.iterrows():
    dA, dB, cell = row['drugA_name'], row['drugB_name'], row['cell_line']

    ic50_A = ic50_map.get((dA, cell))
    ic50_B = ic50_map.get((dB, cell))

    if ic50_A is None or ic50_B is None:
      continue

    mask = (
        (final_dataset['drugA_name'] == dA)
        & (final_dataset['drugB_name'] == dB)
        & (final_dataset['cell_line'] == cell)
    )
    sub_df = final_dataset[mask].copy()

    sub_df['drugA Conc (µM)'] = sub_df['drugA Conc (µM)'].round(4)
    sub_df['drugB Conc (µM)'] = sub_df['drugB Conc (µM)'].round(4)
    sub_df['inhibition'] = 1 - sub_df['X/X0']

    matrix_inhibition = sub_df.pivot_table(
        index='drugA Conc (µM)',
        columns='drugB Conc (µM)',
        values='inhibition',
        aggfunc='mean',
    )
    matrix_inhibition = matrix_inhibition.fillna(0)

    key = f'{dA}__{dB}__{cell}'
    inhibition_matrices[key] = matrix_inhibition

  return inhibition_matrices

def print_dataset_statistics(final_dataset):
  """Print descriptive statistics of the filtered synergy dataset."""
  total_rows = len(final_dataset)

  # Create a unique combination identifier
  final_dataset['combo_id'] = (
      final_dataset['drugA_name']
      + ' + '
      + final_dataset['drugB_name']
      + ' ('
      + final_dataset['cell_line']
      + ')'
  )
  unique_combos = final_dataset['combo_id'].nunique()
  points_per_combo = final_dataset['combo_id'].value_counts()

  print('=' * 50)
  print(f'Total remaining rows (experimental data points): {total_rows}')
  print(f'Total unique drug combinations (possible matrices): {unique_combos}')
  print(f'Average data points per combination: {points_per_combo.mean():.2f}')
  print('=' * 50)

  # Print details to match the Google Colab output:
  print('\nTop drug combinations with the most data points:')
  print(points_per_combo.head(10))
  print('-' * 50)

  print('\nDistribution of data points across matrices:')
  print(points_per_combo.value_counts().sort_index(ascending=False))
  print('=' * 50)