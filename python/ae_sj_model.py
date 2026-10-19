"""ae_sj_model.py - the fitted AE-SJ model and its scoring function.

The model is fitted by R/model.R, which writes it to
outputs/model/FINAL_AE_SJ_covadj_model.json; this module reads that file, so
Python and R apply exactly the same model. After refitting R/model.R nothing
here needs changing.

Scoring:
    1. raw 64-D AlphaEarth embedding (bands A00..A63)
    2. X_adj = AE_COVARIANCE_ADJUSTMENT @ X          (64x64 . 64x1)
    3. X_unit = X_adj / ||X_adj||2                    (L2 normalise)
    4. raw_score = X_unit . AE_SJ_DIRECTION_COVADJ     (dot product)
    5. AE_SJ = (raw_score - SJ4_RAW_ANCHOR) / (SJ1_RAW_ANCHOR - SJ4_RAW_ANCHOR)

0 = median SJ4 reference (late regrowth succession); 1 = median SJ1
reference (intact). Score is intentionally not clipped to [0, 1].
"""

import json
from pathlib import Path

import numpy as np

MODEL_FILE = Path(__file__).resolve().parent.parent / "outputs" / "model" / "FINAL_AE_SJ_covadj_model.json"

_M = json.loads(MODEL_FILE.read_text())
AE_BANDS = list(_M["bands"])
AE_COVARIANCE_ADJUSTMENT = _M["adjustment_matrix"]
AE_SJ_DIRECTION_COVADJ = _M["direction"]
SJ4_RAW_ANCHOR = float(_M["sj4_raw_anchor"])
SJ1_RAW_ANCHOR = float(_M["sj1_raw_anchor"])

_W = np.array(AE_COVARIANCE_ADJUSTMENT, dtype=np.float64)
_DIRECTION = np.array(AE_SJ_DIRECTION_COVADJ, dtype=np.float64)

assert AE_BANDS == [f"A{i:02d}" for i in range(64)], "Unexpected AlphaEarth band order"
assert _W.shape == (64, 64), f"Adjustment matrix must be 64x64, got {_W.shape}"
assert _DIRECTION.shape == (64,), f"Direction vector must have 64 elements, got {_DIRECTION.shape}"


def score_ae_sj(embeddings: np.ndarray) -> np.ndarray:
    """Score raw 64-D AlphaEarth embeddings to the AE-SJ axis.

    embeddings: array of shape (..., 64), band order matching AE_BANDS
    (A00..A63). Any leading shape is preserved (e.g. (n_pixels, 64) or
    (rows, cols, 64) for a raster stack).

    Returns an array of the same leading shape with the unclipped AE-SJ
    score. NaN in, NaN out (propagates through matmul/normalisation).
    """
    x = np.asarray(embeddings, dtype=np.float64)
    if x.shape[-1] != 64:
        raise ValueError(f"Expected 64 AlphaEarth bands in the last axis, got {x.shape[-1]}")

    # Step 2: 64x64 . 64x1 for every pixel - x @ W.T applies the same
    # matrix to every row (equivalent to W @ x per pixel).
    x_adj = x @ _W.T

    # Step 3: L2 normalise. Pixels with (near-)zero norm become NaN, not 0 -
    # dividing by ~0 would otherwise silently produce a huge, meaningless
    # score instead of a clearly-missing one.
    norm = np.sqrt(np.sum(x_adj ** 2, axis=-1, keepdims=True))
    with np.errstate(invalid="ignore", divide="ignore"):
        x_unit = np.where(norm > 1e-12, x_adj / norm, np.nan)

    # Step 4-5: project onto the SJ direction, then rescale to the SJ4=0/SJ1=1 axis.
    raw_score = np.sum(x_unit * _DIRECTION, axis=-1)
    ae_sj = (raw_score - SJ4_RAW_ANCHOR) / (SJ1_RAW_ANCHOR - SJ4_RAW_ANCHOR)
    return ae_sj
