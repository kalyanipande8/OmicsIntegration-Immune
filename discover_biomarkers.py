"""
discover_biomarkers.py

Real statistical biomarker discovery: for every candidate feature across all
4 modalities, computes a Spearman correlation against the simulated vaccine
response score, applies Benjamini-Hochberg FDR correction across all ~90
candidate features, and ranks significant hits. Results are validated against
the injected ground-truth biomarker set (never used in discovery itself) to
report precision/recall/F1 of the discovery procedure.

Writes:
  output/biomarker_discovery_results.csv  - all features, stats, FDR q-values
  output/discovered_biomarkers.json       - the significant subset, ranked
  Also upserts a `biomarkers` table into the SQLite store for the API layer.

Run with: python discover_biomarkers.py  (after simulate_multiomics_data.py)
"""

import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = Path("data")
OUTPUT_DIR = Path("output")
DB_PATH = DATA_DIR / "omics.db"
FDR_ALPHA = 0.05


def benjamini_hochberg(pvals):
    pvals = np.asarray(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p-value down
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    samples = pd.read_sql("SELECT * FROM samples", conn)
    modality_tables = {
        "scrna": "modality_scrna",
        "spatial": "modality_spatial",
        "flow": "modality_flow",
        "tcr_bcr": "modality_tcr_bcr",
    }

    merged = samples[["sample_id", "response_score"]].copy()
    feature_modality = {}
    for modality, table in modality_tables.items():
        df = pd.read_sql(f"SELECT * FROM {table}", conn)
        feat_cols = [c for c in df.columns if c != "sample_id"]
        for c in feat_cols:
            feature_modality[c] = modality
        merged = merged.merge(df, on="sample_id", how="left")

    feature_cols = [c for c in merged.columns if c not in ("sample_id", "response_score")]

    results = []
    for feat in feature_cols:
        rho, p = stats.spearmanr(merged[feat], merged["response_score"])
        results.append({"feature": feat, "modality": feature_modality[feat], "spearman_r": rho, "p_value": p})

    res_df = pd.DataFrame(results)
    res_df["q_value"] = benjamini_hochberg(res_df["p_value"].values)
    res_df["significant"] = res_df["q_value"] < FDR_ALPHA
    res_df = res_df.sort_values("q_value").reset_index(drop=True)
    res_df.to_csv(OUTPUT_DIR / "biomarker_discovery_results.csv", index=False)

    discovered = res_df[res_df["significant"]].copy()
    discovered = discovered.sort_values("spearman_r", key=lambda s: s.abs(), ascending=False)

    # Validate against injected ground truth
    with open(DATA_DIR / "ground_truth_biomarkers.json") as f:
        ground_truth = set(json.load(f).keys())

    discovered_set = set(discovered["feature"])
    tp = discovered_set & ground_truth
    fp = discovered_set - ground_truth
    fn = ground_truth - discovered_set
    precision = len(tp) / len(discovered_set) if discovered_set else 0.0
    recall = len(tp) / len(ground_truth) if ground_truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    validation = {
        "n_candidate_features": len(feature_cols),
        "n_significant_discovered": len(discovered_set),
        "n_ground_truth": len(ground_truth),
        "true_positives": sorted(tp),
        "false_positives": sorted(fp),
        "false_negatives": sorted(fn),
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
    }

    discovered_records = discovered.to_dict(orient="records")
    for rec in discovered_records:
        rec["in_ground_truth"] = rec["feature"] in ground_truth

    with open(OUTPUT_DIR / "discovered_biomarkers.json", "w") as f:
        json.dump({"biomarkers": discovered_records, "validation": validation}, f, indent=2)

    # Persist to SQLite for the API layer
    biomarkers_df = discovered.copy()
    biomarkers_df["in_ground_truth"] = biomarkers_df["feature"].isin(ground_truth)
    biomarkers_df["rank"] = range(1, len(biomarkers_df) + 1)
    biomarkers_df.to_sql("biomarkers", conn, if_exists="replace", index=False)
    conn.commit()
    conn.close()

    print(f"Tested {len(feature_cols)} candidate features across 4 modalities.")
    print(f"Discovered {len(discovered_set)} significant biomarkers (BH-FDR q<{FDR_ALPHA}).")
    print(f"Validation vs. {len(ground_truth)} injected ground-truth biomarkers: "
          f"precision={precision:.2f}, recall={recall:.2f}, F1={f1:.2f}")
    return res_df, validation


if __name__ == "__main__":
    main()
