"""
run_all.py

Runs the full pipeline end-to-end:
  1. Simulate 4 omics modalities across 50+ vaccine studies -> data/omics.db
  2. Statistical biomarker discovery + ground-truth validation -> output/
  3. Smoke-test the researcher-facing Flask API with a real test client

Run with:  python run_all.py
"""

import simulate_multiomics_data
import discover_biomarkers
from api.app import create_app


def main():
    print("--- 1. Simulating multi-omics vaccine study data ---")
    simulate_multiomics_data.simulate()

    print("\n--- 2. Statistical biomarker discovery ---")
    discover_biomarkers.main()

    print("\n--- 3. Smoke-testing the researcher API ---")
    app = create_app()
    client = app.test_client()

    r = client.get("/api/studies")
    print(f"GET /api/studies -> {r.status_code}, {r.get_json()['count']} studies")

    r = client.get("/api/biomarkers")
    print(f"GET /api/biomarkers -> {r.status_code}, {r.get_json()['count']} biomarkers")

    r = client.get("/api/health")
    print(f"GET /api/health -> {r.status_code}")

    print("\nPipeline complete. Run `pytest -q` for the full test suite, "
          "or `python api/app.py` to serve the API on http://localhost:5050")


if __name__ == "__main__":
    main()
