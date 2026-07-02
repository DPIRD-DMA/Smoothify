# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Changed
- Start-point variants are now merged by a per-point median consensus instead of a geometric union. Douglas-Peucker is start-vertex dependent, so the four rotated variants disagree slightly about where a feature's vertices land; the union superimposed those disagreements, whereas the median resamples the variants to a common point count, brings them into phase via an FFT cross-correlation, and takes the coordinate-wise median — a start-invariant consensus that outvotes a single bad variant instead of averaging toward it. This is also faster than the union it replaces (~1.3x end-to-end at the default 3 iterations, ~1.4x at 5, on `examples/Water.gpkg`). The median sits at the consensus interior rather than the outer envelope, so the merged ring is marginally smaller before area preservation restores it (with `preserve_area=True`); with `preserve_area=False` the small reduction is left as-is.
- Faster smoothing with no meaningful change to output. Start-point variants now begin at evenly spaced arc-length positions, computed by interpolating a single new start vertex, rather than rotating to a vertex index on a fully densified ring. The old path segmentized the whole ring only so that index-based rotation would land on evenly spaced points — but every vertex it added was collinear and immediately stripped by the following Douglas-Peucker. Feeding DP just the original vertices plus one cuts its input roughly 3.4x and the simplify step ~3.5x on `examples/Water.gpkg` (total single-core CPU work ~4.4s → ~2.7s, about 1.25x faster end-to-end at the default 3 iterations). The per-point median phase-alignment and Chaikin corner cutting also replace `numpy.roll` with a slice-based cyclic shift (bit-identical, fewer allocations). Output is unchanged within tolerance (per-polygon area drift ≤ 0.01%, worst concave turn unchanged).

### Fixed
- Eliminated the residual "dimple" fold that appeared where the start-point variants disagreed about a shallow feature near the simplify tolerance. Their union superimposed the disagreement as a double peak with a valley that the area-preservation buffer then sharpened into a fold slipping under the concave-turn seal threshold. The per-point median consensus (see Changed) resolves the variants to one boundary instead of overlaying them, removing the artifact (raster-noise fuzz harness: 4 → 0 residual folds; `examples/Water.gpkg`: 8 → 0).

## [0.3.2] - 2026-06-30

### Changed
- Area preservation is faster with no change to output. The buffer-distance search now seeds from the Steiner area expansion `A(d) ≈ A₀ + L·d + π·d²` — solving that quadratic for the initial offset instead of using a linear estimate — and refines with Newton's method using the measured boundary length as the derivative (the rate of area change of an outward offset equals its perimeter), so it no longer brackets the root before solving. This cuts buffer operations per ring from roughly 5–6 to 1–2. On `examples/Water.gpkg` the full single-core pipeline is ~1.3x faster at the default 5 iterations (5.1s → 3.8s), with the speedup scaling up with iteration count and the polygon share of the workload (down to ~1.0x on hole-dominated inputs). Output is unchanged within the area tolerance (per-polygon symmetric difference ≤ 0.005%, area-preservation accuracy identical). Awkward shapes where Newton stalls (pinch-offs or topology changes near the root) fall back to the previous bracketed Brent's-method search, so robustness is unchanged.

## [0.3.1] - 2026-06-15

### Fixed
- Fixed an `AssertionError` crash (`Resulting geometry must be Polygon or LineString. Got <GeometryCollection>`) when smoothing certain valid thin/staircase polygons only a few `segment_length`s across. A rotated start-point variant could self-intersect after simplify + Chaikin where it rounds a neck about one `segment_length` wide; `make_valid()` then repaired that variant into a polygon *plus* a zero-area line filament at the pinch, and `unary_union` kept the 1-D debris, so the variant union became a `GeometryCollection` the rest of the pipeline could not handle. Zero-area debris is now stripped from each repaired variant before the union, so the result stays a clean, valid `Polygon`.
- Improved the sharp-concave-fold repair (introduced in 0.3.0) for thin features that shed significant area during smoothing. The closing (dilate-erode) seal now operates on the area-preserved result and restores area again afterwards, instead of sealing the area-deficient pre-buffer variant union — which the area-preservation buffer could then re-expand into the same cusp. This reduces residual folds on narrow arms.

## [0.3.0] - 2026-06-12

### Changed
- Single-core smoothing is ~3.5x faster on typical raster-derived data (benchmarked on `examples/Water.gpkg`: 6.8s → 1.9s, including the cost of the new fold-repair check below). Output differences are sub-pixel (bounded by the algorithm's own start-point noise floor); area preservation accuracy is unchanged or slightly better. The main changes:
  - Removed the reversed-direction smoothing variants for polygons: Chaikin corner cutting is direction-invariant on closed rings, so they duplicated the forward variants bit-for-bit and only inflated the variant union (output unchanged).
  - Holes are now subtracted in a single `difference` call against their union instead of one `intersection` + `difference` pair per hole (output unchanged).
  - The tiny merge/dissolve buffer now uses mitre joins, which keep corners as single vertices instead of adding ~8 arc vertices per corner (boundary differences at the millimetre scale of the buffer itself).
  - Pre-union Chaikin smoothing of the start-point variants is capped at 2 iterations; detail beyond that was erased by the post-union simplify anyway, while doubling the vertex count entering the expensive union. The final smoothing pass still runs the full `smooth_iterations`.
  - The area-preservation root finder brackets the root from one side using its linear estimate, caches evaluations, and uses a step tolerance derived from the area tolerance via the perimeter, roughly halving the number of buffer operations.
  - Congruent geometries (translated copies of the same shape, common in raster-derived data) are now smoothed once and the result translated to each occurrence; this also reduces work dispatched to parallel workers.

### Added
- New `merge_holes` option (default `True`): holes that touch or nearly touch (e.g. diagonally adjacent raster cells) are joined before smoothing, so they smooth into one coherent opening instead of separate overlapping shapes leaving a fake land bridge. Pass `merge_holes=False` for the previous per-hole behaviour. Mirrors what `merge_multipolygons` does for shells.
- `examples/merge_holes_examples.ipynb`: worked examples of `merge_holes` and its interplay with `merge_collection` (touching holes, overlapping donuts, holes split across features).
- `benchmarks/bench_water.py`: single-core timing/profiling benchmark with baseline output comparison.

### Fixed
- Fixed sharp concave folds in smoothed output on shapes with features about one `segment_length` wide (e.g. a hole with a one-pixel-wide arm, as produced by `merge_holes` joining a small hole to a larger one). The start-point smoothing variants can disagree about such features, leaving a forked slit in their union that the area-preservation shrink sharpens into a cusp. The final result is now checked for sharp concave turns and, when one is found, recomputed from the variant union sealed with a small closing (dilate-erode at `segment_length / 4`).
- Fixed sharp cusps where a smoothed hole crosses the independently smoothed exterior and gets clipped by the hole subtraction (previously up to a 180-degree fold at the tangential crossing points). When detected, the result is repaired with a small opening + closing (at `segment_length / 4`), which removes hair-thin material needles and seals thin slits without visibly moving the boundary.
- Fixed shapes with long straight edges being massively over-rounded (e.g. a large square with a small `segment_length` collapsed into a circle). The simplify steps strip all collinear vertices from straight edges, and Chaikin's corner cuts scale with segment length, so corner rounding grew with edge length instead of `segment_length`. Geometries are now re-segmentized after each simplify step, capping corner rounding at roughly `segment_length` while still smoothing raster staircase artifacts into curves. Applies to all geometry types. Note: outputs are somewhat denser (more vertices) and smoothing is ~20% slower on typical raster-derived data.

## [0.2.3] - 2026-06-02

### Added
- The package now ships inline type hints and a `py.typed` marker, so type checkers (mypy, Pyright/Pylance) pick up `smoothify()`'s signatures and overloads in downstream code.

### Changed
- Invalid geometries (e.g. self-intersecting polygons) are now returned unchanged with a warning instead of crashing with a cryptic error or silently collapsing to an empty geometry. Behaviour is consistent across single geometries, lists/collections, and GeoDataFrames (and regardless of `merge_collection`). Repair them with shapely's `make_valid()` first if you want them smoothed.

## [0.2.2] - 2026-03-24

### Fixed
- Fixed LineString smoothing: self-intersecting lines are no longer run through the polygon-only `make_valid`/`unary_union` step, which could split them into a MultiLineString at crossing points
- Preserve coordinate dimensionality during Chaikin corner cutting instead of forcing 2D output, so smoothed geometries keep their original dimensions

## [0.2.1] - 2026-02-25

### Fixed
- Fixed `TopologyException` crash on thin/elongated polygons by validating smoothed variants before union
- Fixed crash when hole subtraction splits a polygon into a MultiPolygon (e.g. tiny holes relative to segment length)

## [0.2.0] - 2026-02-25

### Fixed
- `smooth_iterations=0` now returns the original input unchanged instead of running the geometry through segmentize/simplify pipeline without smoothing

## [0.1.0] - 2025-11-25

### Added
- Initial public release
- Core smoothing functionality using Chaikin's corner-cutting algorithm
- Support for all Shapely geometry types (Polygon, LineString, MultiPolygon, etc.)
- Automatic segment length detection
- Parallel processing support
- Area preservation for polygons