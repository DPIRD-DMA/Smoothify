"""Output-quality sweep over real-world data across a range of settings.

Runs the full Water.gpkg example dataset through smoothify with a range of
parameter combinations and automatically screens every output for defect
classes seen in the wild:

- invalid or empty geometries
- sharp concave folds/cusps (variant-union slits, clipped-hole cusps) — a
  correctly smoothed boundary is gently curved on the concave side
- gross area loss/gain when area preservation is enabled
- dropped features when no merging is requested

Thresholds are calibrated against this dataset: worst observed concave turn
across all configs is ~90 degrees (a sub-pixel dimple at the tip of a
one-pixel-wide river arm); genuine fold artifacts measure 140-180 degrees.
"""

from pathlib import Path

import geopandas as gpd
import pytest

from smoothify import smoothify

from .test_merge_holes import max_concave_turn

WATER = Path(__file__).parent.parent / "examples" / "Water.gpkg"

CONFIGS = {
    # single knobs off defaults
    "defaults": {},
    "fine_segment": {"segment_length": 4.5},
    "coarse_segment": {"segment_length": 18.0},
    "very_coarse_segment": {"segment_length": 45.0},
    "one_iteration": {"smooth_iterations": 1},
    "five_iterations": {"smooth_iterations": 5},
    "no_area_preservation": {"preserve_area": False},
    "no_merge_collection": {"merge_collection": False},
    "no_merge_holes": {"merge_holes": False},
    "tight_area_tolerance": {"area_tolerance": 0.001},
    "merge_field": {"merge_field": "DN"},
    # interactions: artifacts have historically surfaced where knobs combine
    # (the clipped-hole cusp only showed with merge_collection=False)
    "coarse_no_merge": {"segment_length": 18.0, "merge_collection": False},
    "separate_everything": {"merge_collection": False, "merge_holes": False},
    "minimal_smoothing": {"smooth_iterations": 1, "preserve_area": False},
}

# Sharpest concave turn tolerated anywhere in any output ring. Fold-class
# artifacts (the bugs this guards against) measure 140-180 degrees.
MAX_CONCAVE_TURN = 100.0
# Net dataset area drift. With preserve_area each polygon is restored to
# ~0.01%, and merging/hole-joining moves a little more by design.
MAX_TOTAL_AREA_DRIFT = 0.01
# Without area preservation the dataset legitimately shrinks a few percent:
# Chaikin cuts corners inward and the median variant-merge sits at the
# consensus interior rather than the outer envelope. Only guard gross drift.
MAX_TOTAL_AREA_DRIFT_NO_PRESERVE = 0.08


@pytest.fixture(scope="module")
def water_gdf():
    return gpd.read_file(WATER)


@pytest.mark.slow
@pytest.mark.parametrize("name", CONFIGS)
def test_water_outputs_are_clean(water_gdf, name):
    kwargs = CONFIGS[name]

    result = smoothify(water_gdf, num_cores=1, **kwargs)

    assert len(result) > 0
    if kwargs.get("merge_collection") is False:
        # Without merging, every input feature must come back
        assert len(result) == len(water_gdf)

    assert result.geometry.notna().all(), "null geometry in output"
    assert result.geometry.is_valid.all(), "invalid geometry in output"
    assert (~result.geometry.is_empty).all(), "empty geometry in output"

    worst = max(max_concave_turn(geom) for geom in result.geometry)
    assert worst < MAX_CONCAVE_TURN, (
        f"sharp concave fold in output: {worst:.1f} degrees (config {name!r})"
    )

    in_area = water_gdf.geometry.area.sum()
    out_area = result.geometry.area.sum()
    drift = abs(out_area - in_area) / in_area
    limit = (
        MAX_TOTAL_AREA_DRIFT
        if kwargs.get("preserve_area", True)
        else MAX_TOTAL_AREA_DRIFT_NO_PRESERVE
    )
    assert drift < limit, f"total area drifted {drift:.2%} (config {name!r})"


@pytest.mark.slow
def test_parallel_output_matches_serial(water_gdf):
    """The parallel path (congruence dedup + worker dispatch) must produce
    exactly the same geometries as the serial path."""
    serial = smoothify(water_gdf, num_cores=1)
    parallel = smoothify(water_gdf, num_cores=2)

    assert len(parallel) == len(serial)
    assert all(
        a.equals_exact(b, tolerance=0)
        for a, b in zip(serial.geometry, parallel.geometry, strict=True)
    )
