"""Correlation plot for experimental vs predicted ASO scores."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


def _safe_corr(x: pd.Series, y: pd.Series):
    """Return (pearson_r, pearson_p, spearman_rho, spearman_p) or NaNs.

    Handles small n, NaNs, and constant columns gracefully.
    """
    from scipy import stats

    mask = x.notna() & y.notna()
    if mask.sum() < 3:
        return (np.nan, np.nan, np.nan, np.nan)
    xc, yc = x[mask], y[mask]
    if xc.nunique() < 2 or yc.nunique() < 2:
        return (np.nan, np.nan, np.nan, np.nan)
    r, p = stats.pearsonr(xc, yc)
    rho, rp = stats.spearmanr(xc, yc)
    return (float(r), float(p), float(rho), float(rp))


def correlation_plot(
    matched_df: pd.DataFrame,
    score_columns: list[str],
    measured_col: str,
    out_png: Path,
    exon_col: Optional[str] = "Region (Exon)",
) -> dict:
    """Render side-by-side regression plots per exon, one overlay per score column.

    Each score column is min-max normalized per exon to [0, 1] before plotting
    so they can share a common x-axis. Correlations are computed on the raw
    (un-normalized) values.

    Returns nested stats: {exon: {source: (pearson_r, pearson_p, spearman, p)}}
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    if exon_col and exon_col in matched_df.columns:
        exons = sorted(matched_df[exon_col].dropna().unique().tolist())
        panels: list[tuple[str, pd.DataFrame]] = [
            (str(e), matched_df[matched_df[exon_col] == e].copy()) for e in exons
        ]
    else:
        panels = [("all", matched_df.copy())]

    if not panels:
        raise ValueError("No data to plot.")

    fig, axes = plt.subplots(1, len(panels), figsize=(7 * len(panels), 6), squeeze=False)
    axes_flat = axes[0]

    palette = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
    stats_out: dict[str, dict[str, tuple[float, float, float, float]]] = {}

    for ax, (exon, sub) in zip(axes_flat, panels):
        stats_out[exon] = {}
        if len(sub) < 3:
            ax.set_title(f"{exon} (n={len(sub)}, too few)")
            continue

        for col in score_columns:
            if col not in sub.columns:
                continue
            s = sub[col]
            lo, hi = s.min(), s.max()
            denom = (hi - lo) if hi > lo else 1.0
            sub[f"{col}_norm"] = (s - lo) / denom

        for i, col in enumerate(score_columns):
            if col not in sub.columns:
                continue
            r, p, rho, rp = _safe_corr(sub[col], sub[measured_col])
            stats_out[exon][col] = (r, p, rho, rp)
            label = f"{col} (r={r:.3f}, ρ={rho:.3f})" if not np.isnan(r) else f"{col} (n/a)"
            sns.regplot(
                data=sub,
                x=f"{col}_norm",
                y=measured_col,
                ax=ax,
                color=palette[i % len(palette)],
                label=label,
                scatter_kws={"alpha": 0.7, "edgecolors": "k", "linewidths": 0.5},
            )

        ax.set_xlabel("Predicted score (normalized)")
        ax.set_ylabel(measured_col)
        ax.set_title(f"Exon {exon} (n={len(sub)})")
        ax.legend(fontsize=8)

    plt.tight_layout()
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    return stats_out


def print_correlation_table(stats: dict) -> None:
    """Pretty-print the correlation results as a table."""
    rows: list[dict] = []
    for exon, by_source in stats.items():
        for source, (r, p, rho, rp) in by_source.items():
            rows.append(
                {
                    "Exon": exon,
                    "Source": source,
                    "Pearson r": round(r, 4) if not np.isnan(r) else "n/a",
                    "Pearson p": round(p, 4) if not np.isnan(p) else "n/a",
                    "Spearman ρ": round(rho, 4) if not np.isnan(rho) else "n/a",
                    "Spearman p": round(rp, 4) if not np.isnan(rp) else "n/a",
                }
            )
    if not rows:
        print("No correlation results.")
        return
    df = pd.DataFrame(rows)
    print("\n=== Correlation summary ===")
    print(df.to_string(index=False))
