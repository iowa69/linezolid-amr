"""Backwards-compatible re-export of the internal MLST module.

The legacy ``linezolid_amr.mlst`` API wrapped Seemann's external ``mlst``
binary. Starting with v0.1.4 the project ships its own PubMLST-backed
implementation, which lives in :mod:`linezolid_amr.internal_mlst`. We keep
``linezolid_amr.mlst`` as a thin alias so existing import paths continue
to work.
"""

from linezolid_amr.internal_mlst import (  # noqa: F401
    MlstResult,
    MlstNoMatch,
    MlstUnsupportedOrganism,
    PUBMLST_DBS,
    blastn_available,
    fetch_pubmlst_scheme,
    fetch_pubmlst_schemes,
    infer_organism,
    run_internal_mlst as run_mlst,
    scheme_available,
    all_schemes_available,
)

# Legacy mapping kept for any external code still importing it
MLST_SCHEME_TO_ORGANISM = {info["name"]: org for org, info in PUBMLST_DBS.items()}


def mlst_available() -> bool:
    """True when both blastn is on PATH and all PubMLST schemes are present.

    The legacy semantics ("is Seemann's mlst binary on PATH") no longer apply
    because we now ship our own MLST. We return True only when the in-house
    implementation can actually run.
    """
    return blastn_available() and all_schemes_available()
