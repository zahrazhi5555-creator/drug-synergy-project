# Drug Synergy & Asymmetric Game Theory Modeling

A computational bioinformatics framework for evaluating optimal drug combination doses and synergy in cancer cell lines using **Asymmetric Game Theory**, **Nash Equilibrium analysis**, and **Dose-Penalty Optimization**.

---

##  Project Overview
Combination therapy is a cornerstone of modern cancer treatment, yet identifying optimal dose ratios remains a major challenge. This repository implements an advanced game-theoretic model where two interacting drugs act as strategic players aiming to maximize therapeutic inhibition while being penalized for excessive dosing relative to their individual half-maximal inhibitory concentration ($\text{IC}_{50}$) values.

---

##  Dataset & Scope
* **Total Records:** 1,313 experimental records
* **Drugs Analyzed:** 17 distinct pharmacological agents (including Lapatinib, BEZ-235, Etoposide, Dasatinib, 5-FU, AZD1775, geldanamycin, Temozolomide, etc.)
* **Cancer Cell Lines:** 6 key lines evaluated, with primary focus on **A375** (melanoma) and **HT29** (colorectal adenocarcinoma).

---

##  Key Methodological Components
1. **Empirical Reference Concentration ($\text{IC}_{50}$):** Integrated database screening and empirical reference mapping to normalize dose responses.
2. **Asymmetric Payoff Matrix Formulation:** 
   $$\text{Payoff}_A = \text{Inhibition} - \beta \left(\frac{\text{Dose}_A}{\text{IC}_{50, A}}\right)$$
   $$\text{Payoff}_B = \text{Inhibition} - \beta \left(\frac{\text{Dose}_B}{\text{IC}_{50, B}}\right)$$
   where $\beta$ represents the penalty factor for high-dose toxicity.
3. **Equilibrium Extraction:** Identifies **Pure Nash Equilibria** or falls back to **Pareto Optimal** strategies (maximizing total joint payoff) when strict pure equilibria are absent.
4. **Robustness & Stability Analysis:** Evaluates sensitivity across a wide range of penalty factors ($\beta \in [0.0, 2.0]$) to ensure model stability.

---

##  Selected Key Drug Combinations (Featured in Manuscript)
The core analysis focuses on 6 biologically relevant combinations selected based on literature and pharmacological significance:
1. **Lapatinib + BEZ-235** on **HT29** (Matrix size: $2 \times 3$)
2. **Etoposide + Dasatinib** on **HT29** (Matrix size: $3 \times 4$)
3. **5-FU + BEZ-235** on **HT29** (Matrix size: $4 \times 4$)
4. **Dasatinib + BEZ-235** on **A375** (Matrix size: $2 \times 2$)
5. **AZD1775 + geldanamycin** on **A375** (Matrix size: $3 \times 3$)
6. **BEZ-235 + Temozolomide** on **A375** (Matrix size: $2 \times 4$, transposed)

---

##  Repository Structure
```text
drug-synergy-project/
│
├── data/                      # Raw and processed experimental datasets
├── results/                   # Generated summary tables, Excel files, and stability plots
├── src/                       # Source code modules
│   ├── preprocessing.py       # Data loading and IC50 mapping
│   ├── game_theory.py         # Core game theory engine & payoff calculations
│   └── selected_games_analysis.py # Script for analyzing the 6 featured manuscript games
├── main.py                    # Main pipeline orchestrator
├── requirements.txt           # Python dependencies
└── README.md                  # Project documentation

Getting Started:
__1. Prerequisites__
Ensure you have Python 3.8+ installed along with the required libraries:
pip install -r requirements.txt
__2. Running the Analysis__
To run the primary workflow and generate the game-theoretic evaluations and summary tables in your terminal, execute:
python main.py
__3.Key Outputs__
**Payoff Matrices:** Detailed matrix evaluations displaying individual utilities $(U_A, U_B)$ and marking optimal strategies ($\text{★}$).
**Summary Tables:** Exported CSV and Excel files (Selected_6_Games_Summary.xlsx) capturing optimal doses, individual payoffs, total joint payoff, and equilibrium types.
**Robustness Plots:** Sensitivity curves illustrating optimal dose stability across varying penalty factors ($\beta$).
**License & Citation**
This project is developed for academic research. If you use or adapt this code for your research, please cite the corresponding publication.