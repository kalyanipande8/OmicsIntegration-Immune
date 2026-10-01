"""
simulate_multiomics_data.py

Simulates 4 distinct immune-profiling modalities across 50+ vaccine studies
and consolidates everything into one unified SQLite store (data/omics.db).

Modalities (all per-subject, per-study):
  1. scRNA-seq       - pseudobulk expression of 40 immune genes
  2. Spatial (CODEX/Xenium-style) - 15 tissue niche / cell-colocalization scores
  3. Flow cytometry  - 20 immune cell subset frequencies (%)
  4. TCR/BCR repertoire - 15 clonality / diversity / clonotype-frequency features

A ground-truth set of exactly 18 features (spread across all 4 modalities) is
wired, with known linear weights, into a simulated continuous vaccine-response
score (e.g. a neutralizing-antibody titer proxy). All other features are pure
noise, uncorrelated with response. This ground truth is saved so the
biomarker-discovery step can be validated against it rather than graded on
its own homework.

Run with: python simulate_multiomics_data.py
"""

import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

RNG_SEED = 42
N_STUDIES = 55
SUBJECTS_PER_STUDY = (8, 16)  # inclusive range, randomized per study

DATA_DIR = Path("data")
DB_PATH = DATA_DIR / "omics.db"

# ---------------------------------------------------------------------------
# Feature definitions
# ---------------------------------------------------------------------------

SCRNA_GENES = [f"GENE_{g}" for g in [
    "IFNG", "TNF", "IL2", "IL6", "IL10", "IL21", "TBX21", "GATA3", "FOXP3",
    "GZMB", "PRF1", "CD38", "CD69", "CXCR5", "PDCD1", "ICOS", "BCL6", "AICDA",
    "MKI67", "CCR7", "SELL", "IL7R", "KLRG1", "TOX", "TCF7", "CD27", "CD28",
    "HLA_DR", "STAT1", "STAT3", "IRF4", "BATF", "XBP1", "PRDM1", "MS4A1",
    "CD19", "IGHG1", "IGHM", "JCHAIN", "TNFRSF17",
]]

SPATIAL_FEATURES = [
    "GC_density_LN", "Tfh_GC_colocalization", "Bcell_follicle_area",
    "Tcell_zone_density", "DC_Tcell_contact_score", "plasma_cell_medullary_cord",
    "niche_diversity_index", "macrophage_subcapsular_density",
    "Tfh_Bcell_distance_um", "germinal_center_count", "stromal_FDC_density",
    "CD8_infiltration_score", "neutrophil_margin_density",
    "vascular_density_LN", "niche_entropy",
]

FLOW_FEATURES = [
    "pct_CD4_T", "pct_CD8_T", "pct_Treg", "pct_Tfh_circulating", "pct_Bmem",
    "pct_plasmablast", "pct_naive_B", "pct_NK", "pct_monocyte_classical",
    "pct_monocyte_nonclassical", "pct_mDC", "pct_pDC", "pct_Th1", "pct_Th2",
    "pct_Th17", "CD4_CD8_ratio", "pct_activated_CD38_HLADR_CD8",
    "pct_ICOS_Tfh", "pct_switched_memory_B", "pct_granzymeB_CD8",
]

TCR_BCR_FEATURES = [
    "TCR_clonality", "BCR_clonality", "TCR_shannon_diversity",
    "BCR_shannon_diversity", "top10_TCR_clone_freq", "top10_BCR_clone_freq",
    "vaccine_specific_clonotype_freq", "somatic_hypermutation_rate",
    "class_switch_fraction_IgG", "BCR_clone_expansion_index",
    "TCR_public_clonotype_count", "CDR3_length_mean_BCR",
    "convergent_recombination_score", "TCR_Vgene_skew", "BCR_Vgene_skew",
]

ALL_FEATURES = {
    "scrna": SCRNA_GENES,
    "spatial": SPATIAL_FEATURES,
    "flow": FLOW_FEATURES,
    "tcr_bcr": TCR_BCR_FEATURES,
}

# ---------------------------------------------------------------------------
# Ground truth: exactly 18 predictive biomarkers, spread across modalities,
# each with a signed weight used to construct the response score.
# ---------------------------------------------------------------------------

GROUND_TRUTH_BIOMARKERS = {
    # scRNA-seq (5)
    "GENE_IFNG": 1.4, "GENE_TBX21": 1.1, "GENE_CXCR5": 0.9,
    "GENE_BCL6": 0.8, "GENE_MKI67": 0.7,
    # spatial (4)
    "GC_density_LN": 1.3, "Tfh_GC_colocalization": 1.2,
    "germinal_center_count": 1.0, "stromal_FDC_density": 0.6,
    # flow (5)
    "pct_Tfh_circulating": 1.3, "pct_plasmablast": 1.5,
    "pct_activated_CD38_HLADR_CD8": 0.9, "pct_ICOS_Tfh": 1.0,
    "pct_switched_memory_B": 0.8,
    # TCR/BCR (4)
    "vaccine_specific_clonotype_freq": 1.6, "BCR_clonality": 1.0,
    "somatic_hypermutation_rate": 0.9, "class_switch_fraction_IgG": 0.7,
}
assert len(GROUND_TRUTH_BIOMARKERS) == 18

VACCINE_TYPES = ["mRNA", "protein-subunit", "viral-vector", "inactivated", "live-attenuated"]
PLATFORMS_BY_MODALITY_NOTE = "10x Genomics (scRNA), CODEX/Xenium (spatial), Aurora spectral flow, Adaptive immunoSEQ (TCR/BCR)"


def simulate():
    rng = np.random.default_rng(RNG_SEED)
    DATA_DIR.mkdir(exist_ok=True)

    studies = []
    samples = []
    scrna_rows, spatial_rows, flow_rows, tcrbcr_rows = [], [], [], []

    sample_counter = 0
    for s in range(1, N_STUDIES + 1):
        study_id = f"VAX-{s:03d}"
        vaccine_type = rng.choice(VACCINE_TYPES)
        n_subjects = int(rng.integers(SUBJECTS_PER_STUDY[0], SUBJECTS_PER_STUDY[1] + 1))
        studies.append({
            "study_id": study_id,
            "vaccine_type": vaccine_type,
            "n_subjects": n_subjects,
            "phase": rng.choice(["Phase1", "Phase1/2", "Phase2"]),
            "year": int(rng.integers(2016, 2026)),
        })

        # Study-level batch effect so modalities aren't trivially identical
        # across studies (adds realism without touching ground truth signal).
        batch_shift = rng.normal(0, 0.15)

        for subj in range(n_subjects):
            sample_counter += 1
            sample_id = f"{study_id}-S{subj+1:03d}"

            feature_values = {}  # flat dict of every raw feature for this sample

            scrna_vals = {}
            for g in SCRNA_GENES:
                base = rng.normal(5 + batch_shift, 1.0)
                scrna_vals[g] = base
                feature_values[g] = base

            spatial_vals = {}
            for f in SPATIAL_FEATURES:
                base = rng.normal(0.5 + batch_shift * 0.1, 0.15)
                spatial_vals[f] = base
                feature_values[f] = base

            flow_vals = {}
            for f in FLOW_FEATURES:
                base = abs(rng.normal(10 + batch_shift, 5))
                flow_vals[f] = base
                feature_values[f] = base

            tcrbcr_vals = {}
            for f in TCR_BCR_FEATURES:
                base = abs(rng.normal(0.3 + batch_shift * 0.05, 0.1))
                tcrbcr_vals[f] = base
                feature_values[f] = base

            # --- construct ground-truth-driven vaccine response score ---
            z = 0.0
            for feat, weight in GROUND_TRUTH_BIOMARKERS.items():
                val = feature_values[feat]
                mu, sigma = _feature_mu_sigma(feat, rng)
                z += weight * (val - mu) / sigma
            noise = rng.normal(0, 1.8)
            response_score = z + noise
            # squash to a 0-100 "neutralizing titer percentile" proxy
            response_pct = 100 / (1 + np.exp(-response_score / 3))
            responder = int(response_pct >= 50)

            samples.append({
                "sample_id": sample_id,
                "study_id": study_id,
                "subject_index": subj + 1,
                "response_score": round(float(response_pct), 3),
                "responder": responder,
            })

            scrna_rows.append({"sample_id": sample_id, **{k: round(v, 4) for k, v in scrna_vals.items()}})
            spatial_rows.append({"sample_id": sample_id, **{k: round(v, 4) for k, v in spatial_vals.items()}})
            flow_rows.append({"sample_id": sample_id, **{k: round(v, 4) for k, v in flow_vals.items()}})
            tcrbcr_rows.append({"sample_id": sample_id, **{k: round(v, 4) for k, v in tcrbcr_vals.items()}})

    studies_df = pd.DataFrame(studies)
    samples_df = pd.DataFrame(samples)
    scrna_df = pd.DataFrame(scrna_rows)
    spatial_df = pd.DataFrame(spatial_rows)
    flow_df = pd.DataFrame(flow_rows)
    tcrbcr_df = pd.DataFrame(tcrbcr_rows)

    # ------------------------------------------------------------------
    # Consolidate into unified SQLite store
    # ------------------------------------------------------------------
    if DB_PATH.exists():
        DB_PATH.unlink()
    conn = sqlite3.connect(DB_PATH)
    studies_df.to_sql("studies", conn, index=False)
    samples_df.to_sql("samples", conn, index=False)
    scrna_df.to_sql("modality_scrna", conn, index=False)
    spatial_df.to_sql("modality_spatial", conn, index=False)
    flow_df.to_sql("modality_flow", conn, index=False)
    tcrbcr_df.to_sql("modality_tcr_bcr", conn, index=False)

    conn.execute("CREATE INDEX idx_samples_study ON samples(study_id)")
    conn.execute("CREATE INDEX idx_scrna_sample ON modality_scrna(sample_id)")
    conn.execute("CREATE INDEX idx_spatial_sample ON modality_spatial(sample_id)")
    conn.execute("CREATE INDEX idx_flow_sample ON modality_flow(sample_id)")
    conn.execute("CREATE INDEX idx_tcrbcr_sample ON modality_tcr_bcr(sample_id)")
    conn.commit()
    conn.close()

    # Save ground truth (used only for validation, never for discovery)
    with open(DATA_DIR / "ground_truth_biomarkers.json", "w") as f:
        json.dump(GROUND_TRUTH_BIOMARKERS, f, indent=2)

    print(f"Simulated {len(studies_df)} studies, {len(samples_df)} samples.")
    print(f"Modalities: scRNA ({len(SCRNA_GENES)} genes), spatial ({len(SPATIAL_FEATURES)} features), "
          f"flow ({len(FLOW_FEATURES)} features), TCR/BCR ({len(TCR_BCR_FEATURES)} features).")
    print(f"Unified SQLite store written to {DB_PATH}")
    return studies_df, samples_df


def _feature_mu_sigma(feat, rng):
    """Approximate per-feature center/scale used only to standardize the
    ground-truth contribution to the response score at generation time."""
    if feat in SCRNA_GENES:
        return 5.0, 1.0
    if feat in SPATIAL_FEATURES:
        return 0.5, 0.15
    if feat in FLOW_FEATURES:
        return 10.0, 5.0
    return 0.3, 0.1


if __name__ == "__main__":
    simulate()
