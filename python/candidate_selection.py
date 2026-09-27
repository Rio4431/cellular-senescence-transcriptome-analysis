"""Candidate selection logic used in the portfolio analysis."""
import numpy as np
import pandas as pd


def add_recovery_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Add Treatment/Normal ratios and log2-distance-based recovery metrics."""
    out = df.copy()
    out["Ratio_Treatment_Normal"] = out["Treatment Expression value"] / out["Normal Expression value"]
    out["Log2Ratio_Treatment_Normal"] = np.log2(out["Ratio_Treatment_Normal"])
    out["Distance_Senescence"] = out["Log2Ratio_A"].abs()
    out["Distance_Treatment"] = out["Log2Ratio_Treatment_Normal"].abs()
    out["Closer_to_Normal"] = out["Distance_Treatment"] < out["Distance_Senescence"]
    return out


def select_candidates(df: pd.DataFrame) -> dict[str, object]:
    """Reproduce the final 718 -> 330 -> 181 -> 51 candidate-selection workflow."""
    sen_changed = df.loc[(df["Ratio_A"] >= 2) | (df["Ratio_A"] <= 0.5)].copy()
    sen_changed = sen_changed.dropna(subset=["Gene Symbol"]).copy()
    sen_changed = add_recovery_metrics(sen_changed)
    sen_changed["Senescence_direction"] = np.where(sen_changed["Ratio_A"] >= 2, "Up", "Down")

    closer_status = sen_changed.groupby("Gene Symbol")["Closer_to_Normal"].agg(["nunique", "first"])
    closer_genes = set(closer_status.index[(closer_status["nunique"] == 1) & closer_status["first"]])
    inconsistent_closer = set(closer_status.index[closer_status["nunique"] > 1])

    closer = sen_changed.loc[sen_changed["Gene Symbol"].isin(closer_genes)].copy()
    closer["Within_2fold_Normal"] = closer["Ratio_Treatment_Normal"].between(0.5, 2.0, inclusive="both")
    within_status = closer.groupby("Gene Symbol")["Within_2fold_Normal"].agg(["nunique", "first"])
    final_genes = set(within_status.index[(within_status["nunique"] == 1) & within_status["first"]])
    inconsistent_within = set(within_status.index[within_status["nunique"] > 1])

    final = sen_changed.loc[sen_changed["Gene Symbol"].isin(final_genes)].copy()
    final["Recovery_percent"] = (
        (final["Distance_Senescence"] - final["Distance_Treatment"])
        / final["Distance_Senescence"] * 100
    )
    final["Treatment_2fold_change"] = (final["Ratio_B"] >= 2) | (final["Ratio_B"] <= 0.5)

    treatment_status = final.groupby("Gene Symbol")["Treatment_2fold_change"].agg(["nunique", "all"])
    strong_genes = set(treatment_status.index[(treatment_status["nunique"] == 1) & treatment_status["all"]])
    strong = final.loc[final["Gene Symbol"].isin(strong_genes)].copy()

    keep = [
        "Gene Symbol", "Probe Set ID", "Normal Expression value",
        "Senescence Expression value_1", "Treatment Expression value",
        "Senescence_direction", "Ratio_A", "Ratio_B", "Ratio_Treatment_Normal",
        "Recovery_percent", "Treatment_2fold_change", "GO Biological Process Term",
        "GO Cellular Component Term", "GO Molecular Function Term",
    ]
    final_table = final[keep].sort_values(["Senescence_direction", "Recovery_percent"], ascending=[True, False]).reset_index(drop=True)
    strong_table = strong[keep].sort_values("Recovery_percent", ascending=False).reset_index(drop=True)

    up_genes = set(sen_changed.loc[sen_changed["Ratio_A"] >= 2, "Gene Symbol"])
    down_genes = set(sen_changed.loc[sen_changed["Ratio_A"] <= 0.5, "Gene Symbol"])

    summary = {
        "senescence_changed": sen_changed["Gene Symbol"].nunique(),
        "senescence_up": len(up_genes),
        "senescence_down": len(down_genes),
        "closer_to_normal": len(closer_genes),
        "final_candidates": len(final_genes),
        "final_up": len(final_genes & up_genes),
        "final_down": len(final_genes & down_genes),
        "strong_candidates": len(strong_genes),
        "closer_probe_inconsistent": sorted(inconsistent_closer),
        "within_probe_inconsistent": sorted(inconsistent_within),
    }
    return {
        "senescence_changed_probes": sen_changed,
        "closer_probes": closer,
        "final_candidates": final_table,
        "strong_candidates": strong_table,
        "summary": summary,
    }
