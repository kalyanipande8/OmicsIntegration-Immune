"""
api/app.py

Researcher-facing Flask REST API over the consolidated multi-omics vaccine
study SQLite store. This is the distinguishing deliverable of this project:
a reusable access layer other researchers/tools can query, rather than a
one-off analysis notebook.

Endpoints:
  GET /api/studies                         - list all studies (filterable by vaccine_type)
  GET /api/studies/<study_id>               - single study + its samples
  GET /api/samples                          - list samples (filterable by study_id, responder)
  GET /api/samples/<sample_id>              - single sample's metadata
  GET /api/modality/<modality>/<sample_id>  - per-modality raw feature data for a sample
  GET /api/biomarkers                       - ranked, discovered biomarkers (filterable by modality)
  GET /api/biomarkers/<feature>             - single biomarker's stats
  GET /api/query/cross_modality/<sample_id> - all 4 modalities + response, joined, for one sample

Run directly with `python api/app.py` for a dev server, or import
`create_app()` for tests (used by tests/test_pipeline.py via Flask's test
client - no live server required for CI).
"""

import sqlite3
from pathlib import Path

from flask import Flask, g, jsonify, request

DB_PATH = Path(__file__).resolve().parent.parent / "data" / "omics.db"
MODALITY_TABLES = {
    "scrna": "modality_scrna",
    "spatial": "modality_spatial",
    "flow": "modality_flow",
    "tcr_bcr": "modality_tcr_bcr",
}


def get_db():
    if "db" not in g:
        conn = sqlite3.connect(DB_PATH)
        conn.row_factory = sqlite3.Row
        g.db = conn
    return g.db


def create_app(db_path=None):
    app = Flask(__name__)
    if db_path:
        app.config["DB_PATH"] = str(db_path)

    @app.teardown_appcontext
    def close_db(exception=None):
        db = g.pop("db", None)
        if db is not None:
            db.close()

    def db_path_for_request():
        return app.config.get("DB_PATH", str(DB_PATH))

    def query(sql, params=()):
        if "db" not in g:
            conn = sqlite3.connect(db_path_for_request())
            conn.row_factory = sqlite3.Row
            g.db = conn
        cur = g.db.execute(sql, params)
        return [dict(row) for row in cur.fetchall()]

    @app.get("/api/studies")
    def list_studies():
        vaccine_type = request.args.get("vaccine_type")
        if vaccine_type:
            rows = query("SELECT * FROM studies WHERE vaccine_type = ?", (vaccine_type,))
        else:
            rows = query("SELECT * FROM studies")
        return jsonify({"count": len(rows), "studies": rows})

    @app.get("/api/studies/<study_id>")
    def get_study(study_id):
        study = query("SELECT * FROM studies WHERE study_id = ?", (study_id,))
        if not study:
            return jsonify({"error": "study not found"}), 404
        samples = query("SELECT * FROM samples WHERE study_id = ?", (study_id,))
        return jsonify({"study": study[0], "samples": samples})

    @app.get("/api/samples")
    def list_samples():
        study_id = request.args.get("study_id")
        responder = request.args.get("responder")
        sql = "SELECT * FROM samples WHERE 1=1"
        params = []
        if study_id:
            sql += " AND study_id = ?"
            params.append(study_id)
        if responder is not None:
            sql += " AND responder = ?"
            params.append(int(responder))
        rows = query(sql, tuple(params))
        return jsonify({"count": len(rows), "samples": rows})

    @app.get("/api/samples/<sample_id>")
    def get_sample(sample_id):
        rows = query("SELECT * FROM samples WHERE sample_id = ?", (sample_id,))
        if not rows:
            return jsonify({"error": "sample not found"}), 404
        return jsonify(rows[0])

    @app.get("/api/modality/<modality>/<sample_id>")
    def get_modality_data(modality, sample_id):
        if modality not in MODALITY_TABLES:
            return jsonify({"error": f"unknown modality '{modality}'",
                             "valid_modalities": list(MODALITY_TABLES)}), 400
        table = MODALITY_TABLES[modality]
        rows = query(f"SELECT * FROM {table} WHERE sample_id = ?", (sample_id,))
        if not rows:
            return jsonify({"error": "sample not found in this modality"}), 404
        return jsonify({"modality": modality, "sample_id": sample_id, "data": rows[0]})

    @app.get("/api/biomarkers")
    def list_biomarkers():
        modality = request.args.get("modality")
        limit = request.args.get("limit", type=int)
        sql = "SELECT * FROM biomarkers WHERE 1=1"
        params = []
        if modality:
            sql += " AND modality = ?"
            params.append(modality)
        sql += " ORDER BY rank ASC"
        if limit:
            sql += " LIMIT ?"
            params.append(limit)
        rows = query(sql, tuple(params))
        return jsonify({"count": len(rows), "biomarkers": rows})

    @app.get("/api/biomarkers/<feature>")
    def get_biomarker(feature):
        rows = query("SELECT * FROM biomarkers WHERE feature = ?", (feature,))
        if not rows:
            return jsonify({"error": "biomarker not found or not significant"}), 404
        return jsonify(rows[0])

    @app.get("/api/query/cross_modality/<sample_id>")
    def cross_modality_query(sample_id):
        sample = query("SELECT * FROM samples WHERE sample_id = ?", (sample_id,))
        if not sample:
            return jsonify({"error": "sample not found"}), 404
        result = {"sample": sample[0]}
        for modality, table in MODALITY_TABLES.items():
            rows = query(f"SELECT * FROM {table} WHERE sample_id = ?", (sample_id,))
            result[modality] = rows[0] if rows else None
        return jsonify(result)

    @app.get("/api/health")
    def health():
        return jsonify({"status": "ok"})

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True, port=5050)
