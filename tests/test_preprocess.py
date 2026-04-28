import json

import polars as pl

from src.preprocess import preprocess_for_cohort_model


def test_preprocess_drops_structural_missing_and_sparse_features(tmp_path):
    input_path = tmp_path / "merged.parquet"
    output_path = tmp_path / "model_ready.parquet"
    manifest_path = tmp_path / "manifest.json"

    pl.DataFrame(
        {
            "patient_id": ["p1", "p2", "p3", "p4"],
            "cohort": ["BRCA", "BRCA", "LUAD", "LUAD"],
            "ENSG_keep": [1.0, 2.0, 3.0, 4.0],
            "cg_structural": [0.1, 0.2, None, None],
            "cg_sparse": [None, 0.2, None, 0.9],
            "1p": [0.0, 0.1, 0.2, 0.3],
        }
    ).write_parquet(input_path)

    preprocess_for_cohort_model(
        input_path=input_path,
        output_path=output_path,
        manifest_path=manifest_path,
        max_missing_rate=0.25,
    )

    out = pl.read_parquet(output_path)
    manifest = json.loads(manifest_path.read_text())

    assert out.columns == [
        "patient_id",
        "cohort",
        "cohort_code",
        "ENSG_keep",
        "1p",
    ]
    assert manifest["label_map"] == {"BRCA": 0, "LUAD": 1}
    assert manifest["feature_columns"] == ["ENSG_keep", "1p"]
    assert manifest["dropped"]["structural_missing_in_any_cohort"] == ["cg_structural"]
    assert manifest["dropped"]["over_missingness_threshold"] == ["cg_sparse"]
