import os
from src.game_theory import (
    solve_common_payoff_games,
    solve_conflicting_payoff_games,
)
from src.preprocessing import (
    build_inhibition_matrices,
    calculate_ic50_map,
    load_datasets,
    print_dataset_statistics,
)


def main():
    print('--- Starting Drug Synergy & Game Theory Pipeline ---')

# Define file paths directly pointing to pre-filtered processed data
synergy_path = os.path.join(
    'data', 'processed', 'filtered_ic50_synergy_data.csv'
)
single_agent_path = os.path.join(
    'data', 'processed', 'filtered_single_agent_data.csv'
)

# Load Pre-filtered Datasets
print('Loading pre-filtered datasets...')
final_dataset, df_single = load_datasets(synergy_path, single_agent_path)

# Print dataset summary statistics
print_dataset_statistics(final_dataset)

# Calculate IC50 Map
print('Calculating IC50 values...')
ic50_map = calculate_ic50_map(df_single)

# Build Inhibition Matrices
print('Building inhibition matrices...')
inhibition_matrices = build_inhibition_matrices(final_dataset, ic50_map)
print(f'Total valid inhibition matrices built: {len(inhibition_matrices)}')

# Ensure results directory exists
os.makedirs('results', exist_ok=True)

# Solve Common Payoff Games
print('Solving common-payoff games...')
df_common = solve_common_payoff_games(inhibition_matrices)
common_out_path = os.path.join(
    'results', 'Common_Payoff_Game_Results.csv'
)
df_common.to_csv(common_out_path, index=False)
print(f'Common-payoff results saved to {common_out_path}')

# Solve Conflicting Payoff Games
print('Solving conflicting-payoff games...')
df_conflicting = solve_conflicting_payoff_games(
    inhibition_matrices, ic50_map, beta=0.05
)
conflicting_out_path = os.path.join(
    'results', 'Conflicting_Payoff_Game_Results.csv'
)
df_conflicting.to_csv(conflicting_out_path, index=False)
print(f'Conflicting-payoff results saved to {conflicting_out_path}')

print('--- Pipeline Completed Successfully ---')


if __name__ == '__main__':
    main()