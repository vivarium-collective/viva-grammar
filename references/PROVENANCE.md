# Reference provenance

## `physicell_hypoxia_cell_rules.csv`
Verbatim copy of the hypoxia model's rule table from the paper's reference
implementation:

- Source: https://github.com/physicell-models/grammar_samples
  `user_projects/hypoxia/config/cell_rules.csv`
- Paper: Johnson, Bergman, …, Fertig, Macklin (2025).
  "Human interpretable grammar encodes multicellular systems biology models
  to democratize virtual cell laboratories." *Cell* 188, 4711–4733.
  https://doi.org/10.1016/j.cell.2025.06.048

8-column PhysiCell rules grammar:
`cell_type, signal, response(increases|decreases), behavior, saturated_value, half_max, hill_power, applies_to_dead`

## Baseline phenotype values (from the hypoxia `PhysiCell_settings.xml`)
Used to seed `CellDefinitions` (the base `b_0` / min `b_m` the CSV does not carry):

| cell type    | cycle entry (1/min) | apoptosis (1/min) | necrosis (1/min) | migration speed (µm/min) | persistence (min) | bias | transform→motile (1/min) |
|--------------|---------------------|-------------------|------------------|--------------------------|-------------------|------|--------------------------|
| tumor        | 7.2e-4              | 5.31667e-5        | 2.8e-3           | 0.0                      | 0                 | 0.0  | 0.001                    |
| motile tumor | 7.2e-4              | 5.31667e-5        | 2.8e-3           | 0.47                     | 15                | 0.18 | 0 (reversion via rule)   |

Oxygen field: D = 1e5 µm²/min, decay 0.1/min, boundary/initial 38 mmHg,
domain [-750,750]² µm, dx = 20 µm, dt_phenotype = 6 min, max_time = 7200 min (5 days).
