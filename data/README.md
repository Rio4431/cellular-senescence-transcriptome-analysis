# Data

The input workbook is private and is **not distributed in this repository**.

To reproduce the analysis locally, place the workbook at:

`data/expression_data.xlsx`

The analysis reads sheet `Sheet1`. The required fields are:

- `Probe Set ID`
- `Normal Expression value`
- `Senescence Expression value_1`
- treatment expression values (normalized internally to `Treatment Expression value`)
- `Ratio_A` and `Log2Ratio_A` for Senescence / Normal
- `Ratio_B` and `Log2Ratio_B` for Treatment / Senescence
- `Gene Symbol`
- `GO Biological Process Term`
- `GO Cellular Component Term`
- `GO Molecular Function Term`

The raw workbook must not be committed. The repository `.gitignore` excludes spreadsheet and tabular files under `data/`.