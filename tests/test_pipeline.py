"""
tests/test_pipeline.py

Runs the full pipeline against fresh simulated data (fixed seed, so
deterministic) and checks:
  - 50+ studies and a sane number of samples are consolidated
  - biomarker discovery finds ~18 significant, validated biomarkers
  - the researcher-facing Flask API works end-to-end via its test client
"""

import json
import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import simulate_multiomics_data as sim
import discover_biomarkers as disc
from api.app import create_app


@pytest.fixture(scope="module", autouse=True)
def built_pipeline():
    sim.simulate()
    disc.main()
    yield


# ---------------------------------------------------------------------------
# Data consolidation
# ---------------------------------------------------------------------------

def test_at_least_50_studies():
    conn = sqlite3.connect(sim.DB_PATH)
    n_studies = conn.execute("SELECT COUNT(*) FROM studies").fetchone()[0]
    conn.close()
    assert n_studies >= 50


def test_samples_present_and_linked_to_studies():
    conn = sqlite3.connect(sim.DB_PATH)
    n_samples = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    orphans = conn.execute(
        "SELECT COUNT(*) FROM samples s LEFT JOIN studies t ON s.study_id = t.study_id "
        "WHERE t.study_id IS NULL"
    ).fetchone()[0]
    conn.close()
    assert n_samples > 300  # 55 studies x 8-16 subjects
    assert orphans == 0


def test_all_four_modalities_populated():
    conn = sqlite3.connect(sim.DB_PATH)
    for table in ["modality_scrna", "modality_spatial", "modality_flow", "modality_tcr_bcr"]:
        n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        assert n > 0, f"{table} is empty"
    conn.close()


# ---------------------------------------------------------------------------
# Biomarker discovery
# ---------------------------------------------------------------------------

def test_approximately_18_biomarkers_discovered():
    with open(ROOT / "output" / "discovered_biomarkers.json") as f:
        result = json.load(f)
    n_discovered = len(result["biomarkers"])
    # allow a realistic band around the injected ground truth of 18
    assert 14 <= n_discovered <= 26


def test_discovery_validated_against_ground_truth():
    with open(ROOT / "output" / "discovered_biomarkers.json") as f:
        result = json.load(f)
    validation = result["validation"]
    assert validation["n_ground_truth"] == 18
    # discovery procedure should recover most/all true signal with reasonable precision
    assert validation["recall"] >= 0.8
    assert validation["precision"] >= 0.6
    assert validation["f1"] >= 0.7


def test_biomarkers_span_multiple_modalities():
    with open(ROOT / "output" / "discovered_biomarkers.json") as f:
        result = json.load(f)
    modalities = {b["modality"] for b in result["biomarkers"]}
    assert len(modalities) >= 3


# ---------------------------------------------------------------------------
# Researcher API (real Flask test client, no live server needed)
# ---------------------------------------------------------------------------

@pytest.fixture
def client():
    app = create_app()
    app.testing = True
    return app.test_client()


def test_api_health(client):
    r = client.get("/api/health")
    assert r.status_code == 200
    assert r.get_json()["status"] == "ok"


def test_api_list_studies(client):
    r = client.get("/api/studies")
    assert r.status_code == 200
    body = r.get_json()
    assert body["count"] >= 50
    assert body["studies"][0]["study_id"].startswith("VAX-")


def test_api_get_single_study(client):
    r = client.get("/api/studies")
    first_study_id = r.get_json()["studies"][0]["study_id"]
    r2 = client.get(f"/api/studies/{first_study_id}")
    assert r2.status_code == 200
    body = r2.get_json()
    assert body["study"]["study_id"] == first_study_id
    assert len(body["samples"]) == body["study"]["n_subjects"]


def test_api_study_not_found(client):
    r = client.get("/api/studies/NOPE-999")
    assert r.status_code == 404


def test_api_list_samples_filtered_by_study(client):
    studies = client.get("/api/studies").get_json()["studies"]
    sid = studies[0]["study_id"]
    r = client.get(f"/api/samples?study_id={sid}")
    assert r.status_code == 200
    samples = r.get_json()["samples"]
    assert len(samples) > 0
    assert all(s["study_id"] == sid for s in samples)


def test_api_per_modality_retrieval(client):
    sample_id = client.get("/api/samples").get_json()["samples"][0]["sample_id"]
    for modality in ["scrna", "spatial", "flow", "tcr_bcr"]:
        r = client.get(f"/api/modality/{modality}/{sample_id}")
        assert r.status_code == 200
        body = r.get_json()
        assert body["sample_id"] == sample_id
        assert len(body["data"]) > 1  # more than just sample_id key


def test_api_unknown_modality_rejected(client):
    sample_id = client.get("/api/samples").get_json()["samples"][0]["sample_id"]
    r = client.get(f"/api/modality/not_a_modality/{sample_id}")
    assert r.status_code == 400


def test_api_biomarker_ranking_endpoint(client):
    r = client.get("/api/biomarkers")
    assert r.status_code == 200
    body = r.get_json()
    assert 14 <= body["count"] <= 26
    ranks = [b["rank"] for b in body["biomarkers"]]
    assert ranks == sorted(ranks)


def test_api_biomarker_lookup_by_feature(client):
    top_feature = client.get("/api/biomarkers").get_json()["biomarkers"][0]["feature"]
    r = client.get(f"/api/biomarkers/{top_feature}")
    assert r.status_code == 200
    assert r.get_json()["feature"] == top_feature


def test_api_biomarker_not_found(client):
    r = client.get("/api/biomarkers/NOT_A_REAL_GENE")
    assert r.status_code == 404


def test_api_cross_modality_query(client):
    sample_id = client.get("/api/samples").get_json()["samples"][0]["sample_id"]
    r = client.get(f"/api/query/cross_modality/{sample_id}")
    assert r.status_code == 200
    body = r.get_json()
    assert body["sample"]["sample_id"] == sample_id
    for modality in ["scrna", "spatial", "flow", "tcr_bcr"]:
        assert body[modality] is not None
