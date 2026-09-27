"""GO Biological Process enrichment and recovery tests."""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from .preprocessing import explode_go_bp


def benjamini_hochberg(p_values) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment preserving the original analysis method."""
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    m = len(ranked)
    adjusted_ranked = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    adjusted_ranked = np.minimum(adjusted_ranked, 1.0)
    adjusted = np.empty(m, dtype=float)
    adjusted[order] = adjusted_ranked
    return adjusted


def prepare_go_data(df: pd.DataFrame, senescence_changed: pd.DataFrame):
    """Prepare gene-level GO BP annotations for background and senescence-changed genes."""
    background = explode_go_bp(df)[["Gene Symbol", "GO_BP"]].drop_duplicates()
    sen = explode_go_bp(senescence_changed, ["Senescence_direction"])
    sen = sen[["Gene Symbol", "Senescence_direction", "GO_BP"]].drop_duplicates()
    return background, sen


def enrichment_by_direction(background: pd.DataFrame, sen_go: pd.DataFrame, direction: str) -> pd.DataFrame:
    """One-sided Fisher enrichment across all GO BP terms in the measured background."""
    background_genes = set(background["Gene Symbol"])
    target_genes = set(sen_go.loc[sen_go["Senescence_direction"] == direction, "Gene Symbol"])
    background_counts = background.groupby("GO_BP")["Gene Symbol"].nunique()
    target_counts = sen_go.loc[sen_go["Senescence_direction"] == direction].groupby("GO_BP")["Gene Symbol"].nunique()
    rows = []
    for term, bg_with in background_counts.items():
        target_with = int(target_counts.get(term, 0))
        target_without = len(target_genes) - target_with
        non_target_with = int(bg_with) - target_with
        non_target_without = len(background_genes) - len(target_genes) - non_target_with
        odds, p = fisher_exact([[target_with, target_without], [non_target_with, non_target_without]], alternative="greater")
        rows.append({
            "GO_BP": term,
            f"{direction}_genes": target_with,
            "Background_genes": int(bg_with),
            "Odds_ratio": odds,
            "P_value": p,
        })
    out = pd.DataFrame(rows).sort_values("P_value").reset_index(drop=True)
    out["FDR"] = benjamini_hochberg(out["P_value"].to_numpy())
    return out


def selective_recovery_test(sen_go: pd.DataFrame, enrichment: pd.DataFrame, final_candidates: pd.DataFrame, direction: str) -> pd.DataFrame:
    """Test whether significant senescence GO terms are selectively represented among recovered candidates."""
    sig_terms = enrichment.loc[enrichment["FDR"] < 0.05, "GO_BP"].tolist()
    annotated = set(sen_go.loc[sen_go["Senescence_direction"] == direction, "Gene Symbol"])
    recovered = annotated & set(final_candidates.loc[final_candidates["Senescence_direction"] == direction, "Gene Symbol"])
    rows = []
    for term in sig_terms:
        term_genes = set(sen_go.loc[(sen_go["Senescence_direction"] == direction) & (sen_go["GO_BP"] == term), "Gene Symbol"])
        a = len(term_genes & recovered)
        b = len(term_genes - recovered)
        nonterm = annotated - term_genes
        c = len(nonterm & recovered)
        d = len(nonterm - recovered)
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"GO_BP": term, "Term_genes": len(term_genes), "Recovered_in_term": a,
                     "Recovery_percent": a / len(term_genes) * 100, "Odds_ratio": odds, "P_value": p})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values("P_value").reset_index(drop=True)
    out["FDR"] = benjamini_hochberg(out["P_value"].to_numpy())
    return out
