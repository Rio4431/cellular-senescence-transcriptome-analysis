# Cellular Senescence Transcriptome Analysis

## Overview

This project is an exploratory microarray analysis of human dermal fibroblasts under three conditions: Normal, Senescence, and Treatment. The workflow identifies genes with large expression changes during senescence and evaluates whether treatment shifts those genes toward the Normal expression state.

The repository is organized as a portfolio project: the notebook presents the analysis story, while reusable data-processing and statistical functions are separated into Python modules.

## Research Question / Objective

1. Which genes change by at least two-fold between Normal and Senescence?
2. Among those genes, which move closer to the Normal expression level after Treatment?
3. Which recovered genes return to within 0.5–2.0-fold of Normal?
4. Which Biological Process terms are over-represented among senescence-associated genes?

## Dataset

The analysis uses private human dermal fibroblast microarray data with three experimental conditions:

- **Normal**: normal dermal fibroblasts
- **Senescence**: experimentally induced senescent fibroblasts
- **Treatment**: senescent fibroblasts after treatment

The raw workbook is not included in this repository. See `data/README.md` for the expected local file location and required fields.

## Methods

- Data validation and Gene Symbol cleanup
- Fold-change filtering: `Ratio >= 2.0` or `Ratio <= 0.5`
- Treatment recovery measured as change in absolute log2 distance from Normal
- Probe-level consistency checks before gene-level classification
- GO Biological Process enrichment using one-sided Fisher's exact tests
- Benjamini-Hochberg false discovery rate correction
- Descriptive visualization of candidate selection and recovery

Because each condition contains a single sample (`n=1`), gene-level inferential differential-expression testing is not performed. Gene candidates are defined using pre-specified fold-change and recovery criteria. GO enrichment tests are gene-set-level analyses and should not be interpreted as gene-level statistical significance.

## Repository Structure

```text
.
├── README.md
├── .gitignore
├── requirements.txt
├── analysis/
│   └── senescence_treatment_analysis.ipynb
├── python/
│   ├── __init__.py
│   ├── preprocessing.py
│   ├── candidate_selection.py
│   └── go_analysis.py
├── data/
│   └── README.md
└── results/
    ├── figures/
    └── tables/
```

## Reproducibility

The analysis was developed for Python/Jupyter workflows in Positron. Dependencies are listed in `requirements.txt`. The private input data are intentionally excluded from version control.

The final workflow preserves the original analysis definitions and checks the expected candidate counts:

`718 -> 330 -> 181 -> 51`

## How to Run

1. Clone the repository.
2. Create and activate a Python environment.
3. Install dependencies with `pip install -r requirements.txt`.
4. Place the private workbook at `data/expression_data.xlsx`.
5. Open `analysis/senescence_treatment_analysis.ipynb` in Positron or Jupyter.
6. Run the notebook from top to bottom.

## Results

- From 24,351 microarray probe sets, 718 genes showed at least a two-fold expression change between Normal and Senescence: 254 Up and 464 Down.
- 330 of these genes moved closer to Normal after Treatment.
- 181 genes were retained as recovery candidates after moving closer to Normal and reaching 0.5–2.0-fold of the Normal expression level.
- Among them, 51 genes were classified as strong Treatment-response candidates because they also changed by at least two-fold between Treatment and Senescence.
- GO Biological Process analysis identified 60 FDR < 0.05 terms among Senescence-Down genes and 4 among Senescence-Up genes.
- No tested significant senescence GO term showed selective recovery at FDR < 0.05 after Treatment.

## Limitations

- Each condition contains one sample (`n=1`), so gene-level statistical significance cannot be evaluated.
- Fold-change thresholds are exploratory selection criteria rather than inferential tests.
- Multiple probes can map to the same Gene Symbol; probe consistency is therefore checked explicitly.
- GO enrichment depends on the annotations represented in the measured microarray background.
- Movement toward the Normal expression level does not establish restoration of biological function or causality.
- The private source data are not distributed, so full external reproduction requires authorized access to the input workbook.
