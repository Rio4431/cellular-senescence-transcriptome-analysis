"""Data loading and preprocessing utilities for the senescence analysis."""
from pathlib import Path
import numpy as np
import pandas as pd

REQUIRED_COLUMNS = [
    "Probe Set ID", "Normal Expression value", "Senescence Expression value_1",
    "Ratio_A", "Log2Ratio_A", "Treatment Expression value", "Ratio_B",
    "Log2Ratio_B", "Gene Symbol", "GO Biological Process Term",
    "GO Cellular Component Term", "GO Molecular Function Term",
]


def load_expression_data(path: str | Path, sheet_name: str = "Sheet1") -> pd.DataFrame:
    """Load the private Excel input and apply minimal text cleanup."""
    df = pd.read_excel(path, sheet_name=sheet_name)
    # The private source workbook uses an internal label for the treatment-expression column.
    # Normalize that sixth column to the public portfolio name without exposing the internal label.
    if "Treatment Expression value" not in df.columns and len(df.columns) > 5:
        df = df.rename(columns={df.columns[5]: "Treatment Expression value"})
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    df = df.copy()
    df["Gene Symbol"] = df["Gene Symbol"].str.strip()
    return df


def ratio_log2_qc(df: pd.DataFrame) -> dict[str, float]:
    """Return maximum absolute differences between ratios and stored log2 ratios."""
    return {
        "Ratio_A_max_abs_diff": float((np.log2(df["Ratio_A"]) - df["Log2Ratio_A"]).abs().max()),
        "Ratio_B_max_abs_diff": float((np.log2(df["Ratio_B"]) - df["Log2Ratio_B"]).abs().max()),
    }


def explode_go_bp(df: pd.DataFrame, extra_columns: list[str] | None = None) -> pd.DataFrame:
    """Convert the delimiter-separated GO Biological Process field to long format."""
    extra_columns = extra_columns or []
    cols = ["Gene Symbol", *extra_columns, "GO Biological Process Term"]
    out = (
        df[cols]
        .dropna(subset=["Gene Symbol", "GO Biological Process Term"])
        .assign(GO_BP=lambda x: x["GO Biological Process Term"].str.split("//"))
        .explode("GO_BP")
    )
    out["GO_BP"] = out["GO_BP"].str.strip()
    return out.loc[out["GO_BP"].ne("")].copy()
