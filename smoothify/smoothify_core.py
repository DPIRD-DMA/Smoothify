import math
from typing import cast

import numpy as np
import numpy.typing as npt
from scipy.optimize import brentq
from shapely import make_valid
from shapely.geometry import (
    GeometryCollection,
    LinearRing,
    LineString,
    MultiPolygon,
    Polygon,
)
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union

# Chaikin's corner cut moves a right-angle vertex inward by ~1/4 of the
# adjacent segment length, so smoothing segments of FACTOR * segment_length
# deviates at most ~segment_length at a right angle — the same deviation
# budget the simplify(segment_length) step already allows — while letting
# shallow facet chains (e.g. simplified curves) blend into smooth curves.
_CHAIKIN_SEGMENT_FACTOR = 4


def _chaikin_corner_cutting(
    geom: Polygon | LineString, num_iterations: int = 1, reverse: bool = False
) -> Polygon | LineString:
    """Apply Chaikin's corner cutting algorithm to smooth a geometry.

    Chaikin's algorithm iteratively replaces each line segment with two new segments
    by cutting corners at 1/4 and 3/4 positions, creating a smooth curve. This
    implementation uses vectorized NumPy operations for performance and handles
    both closed (Polygon) and open (LineString) geometries."""

    is_closed = isinstance(geom, Polygon)
    points = np.array(
        geom.exterior.coords if is_closed else geom.coords, dtype=np.float64
    )

    if is_closed:
        # Remove duplicate endpoint for closed rings
        points = points[:-1]
        endpoints = None
    else:
        # Store endpoints for open linestrings
        endpoints = (points[0], points[-1])

    if reverse:
        points = points[::-1]

    for _ in range(num_iterations):
        # Get point pairs for corner cutting
        if is_closed:
            p0 = points
            p1 = np.roll(points, -1, axis=0)
        else:
            p0 = points[:-1]
            p1 = points[1:]

        # Vectorized smoothing at 1/4 and 3/4 positions
        # Pre-allocate result array for better performance
        n_new_points = len(p0) * 2
        points = np.empty((n_new_points, points.shape[1]), dtype=np.float64)
        points[0::2] = 0.75 * p0 + 0.25 * p1  # q points
        points[1::2] = 0.25 * p0 + 0.75 * p1  # r points

        # Restore endpoints for open linestrings
        if endpoints is not None:
            points = np.vstack([endpoints[0], points, endpoints[1]])

    # Add closing point for polygons
    if is_closed:
        points = np.vstack([points, points[0:1]])

    if reverse:
        points = points[::-1]

    return LineString(points) if isinstance(geom, LineString) else Polygon(points)


def _max_concave_turn_degrees(geom: BaseGeometry) -> float:
    """Sharpest concave (inward) turn angle in degrees across all rings.

    A correctly smoothed polygon is everywhere gently curved on the concave
    side; thin features may legitimately end in sharp convex hairpins, but a
    sharp concave turn means a fold/slit (e.g. from smoothing variants that
    disagreed about a feature near the simplify tolerance)."""

    polys = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
    worst = 0.0
    for poly in polys:
        if not isinstance(poly, Polygon):
            continue
        for ring in [poly.exterior, *poly.interiors]:
            coords = np.asarray(ring.coords)[:-1, :2]
            if len(coords) < 4:
                continue
            vectors = np.diff(np.vstack([coords, coords[:2]]), axis=0)
            cross = vectors[:-1, 0] * vectors[1:, 1] - vectors[:-1, 1] * vectors[1:, 0]
            dot = (vectors[:-1] * vectors[1:]).sum(axis=1)
            # Shoelace sign gives ring orientation; turns with cross sign
            # opposite to it are concave (material on the inside of the bend)
            orientation = np.sign(
                np.dot(coords[:, 0], np.roll(coords[:, 1], -1))
                - np.dot(coords[:, 1], np.roll(coords[:, 0], -1))
            )
            concave = cross * orientation < 0
            if not concave.any():
                continue
            angles = np.degrees(np.arctan2(np.abs(cross[concave]), dot[concave]))
            worst = max(worst, float(angles.max()))
    return worst


def _rotate_ring_coords(ring: LinearRing, shift: float) -> "npt.NDArray[np.float64]":
    """Rotate a linear ring coordinate sequence by a fractional shift.

    Used to create multiple starting point variants of a polygon for smoothing,
    which are later merged to avoid artifacts from a fixed starting vertex."""

    coords = np.array(ring.coords)
    n = len(coords) - 1  # Exclude duplicate closing point
    shift_idx = int(round(shift * n)) % n
    # Use numpy array slicing and concatenation
    rotated = np.vstack(
        [
            coords[shift_idx:],
            coords[1 : shift_idx + 1],
            coords[shift_idx : shift_idx + 1],
        ]
    )
    return rotated


def _rotate_polygon_start(polygon: Polygon, shift: float) -> Polygon:
    """Rotate a polygon's exterior start vertex by a fractional shift.

    Creates a topologically identical polygon with a different starting vertex.
    Used in conjunction with smoothing to avoid artifacts at the start/end point.
    """

    if polygon.is_empty:
        return polygon

    exterior = LinearRing(_rotate_ring_coords(polygon.exterior, shift))
    return Polygon(exterior)


def _generate_starting_point_variants(
    geom: BaseGeometry, n_starting_points: int
) -> list[Polygon | LineString]:
    """Generate evenly rotated variants of a Polygon to avoid smoothing artifacts.

    For Polygons, creates multiple versions with different starting vertices. When
    smoothed and merged, these variants eliminate artifacts that can occur at the
    arbitrary start/end point of a polygon's coordinate ring. LineStrings are
    returned as-is since they have fixed endpoints."""

    if n_starting_points <= 0:
        return []
    if isinstance(geom, LineString):
        return [geom]

    elif isinstance(geom, Polygon):
        if n_starting_points > len(geom.exterior.coords) - 1:
            n_starting_points = len(geom.exterior.coords) - 1

        variants = []

        for i in range(n_starting_points):
            shift = i / n_starting_points
            if shift == 0:
                shifted_geoms = geom
            else:
                shifted_geoms = _rotate_polygon_start(geom, shift)

            variants.append(shifted_geoms)

        return variants
    elif isinstance(geom, MultiPolygon):
        # Invalid input is screened out before smoothing, so a valid Polygon
        # should never become multi-part here. If it does, segmentize/simplify
        # split it unexpectedly — surface that rather than a bare type error.
        raise ValueError(
            "Preprocessing produced a MultiPolygon from a single Polygon "
            f"({len(geom.geoms)} parts); expected it to stay single-part. "
            "This usually means the input geometry was invalid "
            "(self-intersecting). Repair it with shapely's make_valid() first."
        )
    else:
        raise ValueError(
            f"Input geometry must be a Polygon or LineString, got {type(geom)}."
        )


# Most rings land within tolerance in one or two Newton steps; allow a few
# more for awkward shapes before handing off to the bracketed fallback.
_AREA_NEWTON_MAX_STEPS = 6


def _preserve_area_with_buffer(
    polygon: Polygon,
    target_area: float,
    tolerance: float = 1e-6,
) -> Polygon:
    """Restore original polygon area after smoothing via buffering.

    Smoothing slightly changes polygon area; this finds the buffer distance that
    restores the original area to within ``tolerance``.

    The buffered area follows the Steiner expansion ``A(d) ~= A0 + L*d + pi*d^2``
    (``L`` = perimeter), so we seed the search by solving that quadratic instead
    of the cruder linear estimate, then refine with Newton's method. The area
    swept by an outward offset changes at a rate equal to the boundary length
    (the coarea identity), so each buffered candidate yields a near-exact
    derivative for free and no bracketing is needed -- typically one or two
    buffers per ring instead of the half-dozen a bracketed root-find spends.
    Awkward shapes (pinch-offs, topology changes near the root) fall back to the
    robust Brent's-method search."""

    if polygon.is_empty:
        return polygon

    current_area = polygon.area
    if abs(current_area - target_area) <= tolerance:
        return polygon

    perimeter = polygon.length
    if perimeter <= 0:
        return polygon

    # Seed from the Steiner quadratic pi*d^2 + L*d + (A0 - target) = 0, taking
    # the root nearest zero. When the parabola never reaches the target (a deep
    # shrink past its vertex) fall back to the first-order estimate and let
    # Newton walk in.
    area_gap = current_area - target_area
    discriminant = perimeter * perimeter - 4.0 * math.pi * area_gap
    if discriminant >= 0:
        distance = (-perimeter + math.sqrt(discriminant)) / (2.0 * math.pi)
    else:
        distance = -area_gap / perimeter

    best_result: Polygon | None = None
    best_error = float("inf")
    for _ in range(_AREA_NEWTON_MAX_STEPS):
        candidate = polygon.buffer(distance)
        candidate_area = candidate.area
        error = abs(candidate_area - target_area)
        if error < best_error:
            best_result, best_error = candidate, error
        if error <= tolerance:
            return candidate
        # dA/dd == boundary length of the current candidate (coarea identity).
        slope = candidate.length
        if slope <= 0:
            break
        next_distance = distance - (candidate_area - target_area) / slope
        if not math.isfinite(next_distance):
            break
        distance = next_distance

    # Newton stalled before reaching tolerance: defer to the bracketed search,
    # but keep whichever result is actually closer to the target.
    fallback = _preserve_area_brentq(
        polygon, target_area, tolerance, current_area, perimeter
    )
    if fallback is not None and abs(fallback.area - target_area) < best_error:
        return fallback
    return best_result if best_result is not None else polygon


def _preserve_area_brentq(
    polygon: Polygon,
    target_area: float,
    tolerance: float,
    current_area: float,
    perimeter: float,
) -> Polygon | None:
    """Bracketed Brent's-method area restoration (robust fallback).

    Slower than the Newton path (it brackets the root before solving) but does
    not rely on a good initial estimate, so it covers shapes where Newton
    stalls. Returns ``None`` only if no usable buffer could be produced."""

    initial_guess = (target_area - current_area) / perimeter if perimeter > 0 else 0

    # Cache evaluations: brentq re-evaluates the bracket endpoints, and
    # buffer(0) at the known endpoint costs a full (pointless) GEOS buffer.
    _cache: dict[float, float] = {0.0: current_area - target_area}

    def area_delta(distance: float) -> float:
        cached = _cache.get(distance)
        if cached is not None:
            return cached
        result = float(polygon.buffer(distance).area - target_area)
        _cache[distance] = result
        return result

    # Buffered area is monotonic in distance and f(0) = current - target is
    # known for free, so bracket one-sided from the linear estimate instead
    # of buffering both sides of a symmetric interval.
    if initial_guess != 0:
        far = initial_guess * 1.2
        lo, hi = (0.0, far) if far > 0 else (far, 0.0)
    else:
        scale = (polygon.area / 3.1416) ** 0.5
        max_distance = max(0.1, scale)
        lo, hi = -max_distance, max_distance

    f_lo, f_hi = area_delta(lo), area_delta(hi)

    for _ in range(20):
        if f_lo * f_hi < 0:
            break
        lo *= 1.5
        hi *= 1.5
        f_lo, f_hi = area_delta(lo), area_delta(hi)

    try:
        # Use brentq to find optimal buffer distance.
        # xtol is in distance units while tolerance is in area units; since
        # d(area)/d(distance) ~= perimeter, a distance error of
        # tolerance / perimeter corresponds to an area error of ~tolerance.
        # Aim an order of magnitude below that so the verification below
        # rarely needs a second pass.
        xtol = tolerance / perimeter * 0.1 if perimeter > 0 else (hi - lo) * 0.001
        optimal_distance = cast(float, brentq(area_delta, lo, hi, xtol=xtol))
        result = polygon.buffer(optimal_distance)

        # Verify the result meets tolerance, otherwise try to refine
        if abs(result.area - target_area) > tolerance:
            # Try a tighter search if we didn't meet tolerance
            xtol = min(tolerance * 0.001, (hi - lo) * 0.0001)
            optimal_distance = cast(float, brentq(area_delta, lo, hi, xtol=xtol))
            result = polygon.buffer(optimal_distance)

        return result
    except ValueError:
        # Fallback: just return the closest one
        buffered_lo = polygon.buffer(lo)
        buffered_hi = polygon.buffer(hi)
        return min(
            [buffered_lo, buffered_hi],
            key=lambda candidate: abs(candidate.area - target_area),
        )


def _join_adjacent(
    geom: list[BaseGeometry] | BaseGeometry, segment_length: float
) -> BaseGeometry:
    """Merge adjacent geometries using a small buffer-dissolve operation.

    Applies a tiny buffer (segment_length / 1000) to join geometries that are nearly
    touching or share edges, then unions them together into a single geometry or
    collection. This is useful for merging polygons derived from adjacent raster
    cells with the same classification.

    Note: Only applies to Polygon and MultiPolygon geometries. Other geometry types
    (LineString, LinearRing, etc.) are returned unchanged."""

    if isinstance(geom, list):
        geom_combined = unary_union(geom)
    else:
        geom_combined = geom
    poly_types = []
    other_types = []
    if isinstance(geom_combined, GeometryCollection):
        for geom in geom_combined.geoms:
            if isinstance(geom, Polygon | MultiPolygon):
                poly_types.append(geom)
            else:
                other_types.append(geom)
        geom_combined = GeometryCollection(poly_types)
    elif not isinstance(geom_combined, Polygon | MultiPolygon):
        # If the geometry is not a polygon type, return it unchanged
        return geom_combined

    # Mitre join keeps corners as single vertices (round joins add ~8 vertices
    # per corner, inflating the union for no visible gain at this tiny buffer)
    merged_geom = geom_combined.buffer(segment_length / 1000, join_style="mitre")
    merged_geom = unary_union(merged_geom)

    if other_types:
        return GeometryCollection([merged_geom] + other_types)

    return merged_geom


def _polygonal_only(geom: BaseGeometry) -> BaseGeometry:
    """Keep only the polygonal (area-bearing) parts of a geometry.

    A rotated start-point variant can self-intersect after simplify + Chaikin
    where it rounds a neck only ~segment_length wide so the two sides cross.
    make_valid() then repairs it into a Polygon *plus* a zero-area line filament
    at the pinch point. Those 1-D pieces carry no area, but unary_union keeps
    them, so the variant union leaks through as a GeometryCollection (which the
    rest of the pipeline can't handle). Drop them so the union stays polygonal.
    """
    if isinstance(geom, GeometryCollection):
        polys = [g for g in geom.geoms if isinstance(g, (Polygon, MultiPolygon))]
        return unary_union(polys) if polys else Polygon()
    return geom


# Bounds on the per-point resample count used when averaging variants. The count
# scales with the variants' own detail (their max vertex count): too few points
# coarsen a large, detailed ring and make the merged result diverge from the
# per-vertex union, while the floor keeps the phase-alignment accurate on small
# rings (below it, residual folds slip through). The cap bounds the FFT cost on
# pathologically dense rings (some carry tens of thousands of vertices) at a
# resolution well beyond what smoothing needs.
_MIN_RESAMPLE_POINTS = 240
_MAX_RESAMPLE_POINTS = 4000


def _resample_ring(polygon: Polygon, n_points: int) -> "npt.NDArray[np.float64]":
    """Sample a polygon's exterior ring at ``n_points`` equally spaced by arc length."""
    coords = np.asarray(polygon.exterior.coords, dtype=np.float64)
    seg_len = np.sqrt((np.diff(coords, axis=0) ** 2).sum(axis=1))
    cumulative = np.concatenate([[0.0], np.cumsum(seg_len)])
    total = cumulative[-1]
    targets = np.linspace(0.0, total, n_points, endpoint=False)
    x = np.interp(targets, cumulative, coords[:, 0])
    y = np.interp(targets, cumulative, coords[:, 1])
    return np.column_stack([x, y])


def _align_ring(
    sample: "npt.NDArray[np.float64]",
    reference: "npt.NDArray[np.float64]",
    reference_fft: "npt.NDArray[np.complex128]",
) -> "npt.NDArray[np.float64]":
    """Cyclically rotate ``sample`` so its points correspond to ``reference``.

    The variants trace the same shape but start at different vertices, so before
    averaging point-for-point they must be brought into phase. The best cyclic
    shift maximises the cross-correlation, computed here in the complex plane
    (``x + iy``) so a single FFT covers both coordinates; the real part of the
    complex cross-correlation is exactly the summed x/y dot product. The
    reference's transform is precomputed once by the caller. Both roll
    directions are tried to stay agnostic to the FFT sign convention."""
    s = sample - sample.mean(axis=0)
    sample_fft = np.fft.fft(s[:, 0] + 1j * s[:, 1])
    corr = np.fft.ifft(reference_fft * np.conj(sample_fft)).real
    k = int(np.argmax(corr))
    candidates = (np.roll(sample, k, axis=0), np.roll(sample, -k, axis=0))
    return min(candidates, key=lambda c: float(((c - reference) ** 2).sum()))


def _combine_variants(variants: list[BaseGeometry]) -> BaseGeometry:
    """Merge the smoothed start-point variants by a per-point median.

    Douglas-Peucker is start-vertex dependent, so the variants disagree on where
    a feature's vertices land. Taking their geometric *union* superimposes those
    disagreements (a single peak becomes two, with a valley that the later
    area-preservation buffer sharpens into a fold). Averaging instead resolves
    them to consensus, and the median (rather than the mean) outvotes a single
    bad start-point variant rather than being dragged toward it. Each variant is
    resampled to a common point count, brought into phase, and the coordinate-
    wise median is taken; ``make_valid`` cleans up any residual self-touch.

    The median sits at the consensus interior rather than the union's outer
    envelope, so the merged ring is slightly smaller than the union's was; area
    preservation (when enabled) restores it, and when disabled the small
    reduction is acceptable."""
    rings: list[Polygon] = []
    for geom in variants:
        if isinstance(geom, MultiPolygon):
            geom = max(geom.geoms, key=lambda p: p.area)
        if isinstance(geom, Polygon) and not geom.is_empty:
            rings.append(geom)
    if len(rings) < 2:
        return rings[0] if rings else Polygon()

    # Resample every variant to a shared point count so a per-point median is
    # defined. Scale it to the variants' own detail (bounded) so large rings are
    # not coarsened; see _MIN/_MAX_RESAMPLE_POINTS.
    detail = max(len(r.exterior.coords) - 1 for r in rings)
    n_points = int(np.clip(detail, _MIN_RESAMPLE_POINTS, _MAX_RESAMPLE_POINTS))
    samples = [_resample_ring(r, n_points) for r in rings]
    reference = samples[0]
    centered = reference - reference.mean(axis=0)
    reference_fft = np.fft.fft(centered[:, 0] + 1j * centered[:, 1])
    aligned = [reference] + [
        _align_ring(s, reference, reference_fft) for s in samples[1:]
    ]
    combined = np.median(np.array(aligned), axis=0)
    return make_valid(Polygon(np.vstack([combined, combined[:1]])))


def _smoothify_geometry(
    geom: Polygon | LineString,
    segment_length: float,
    smooth_iterations: int = 3,
    preserve_area: bool = True,
    area_tolerance: float = 0.01,
) -> Polygon | LineString:
    """Core smoothing pipeline using Chaikin corner cutting with area preservation.

    This is the main smoothing algorithm that:
    1. Adds intermediate vertices along line segments (segmentize)
    2. Generates multiple rotated variants (for Polygons) to avoid artifacts
    3. Simplifies each variant to remove noise, then re-segmentizes so no
       segment exceeds segment_length * _CHAIKIN_SEGMENT_FACTOR (simplify
       strips collinear vertices, and Chaikin's corner cuts scale with
       segment length — unbounded segments would over-round sharp corners)
    4. Applies Chaikin corner cutting to smooth
    5. Merges the variants by a per-point median (a start-invariant consensus
       that resolves start-point disagreements instead of superimposing them)
    6. Simplifies the merged result and re-segmentizes again, then applies
       a final smoothing pass
    7. Optionally restores original area via buffering (for Polygons)"""

    if geom.geom_type == "Polygon":
        original_area = geom.area
    elif geom.geom_type == "LineString":
        original_area = None
    else:
        raise ValueError("Input geometry must be a Polygon or LineString.")

    geom_segmented = geom.segmentize(segment_length / 2)

    geom_iterations = []

    if isinstance(geom, Polygon):
        starting_point_geoms = _generate_starting_point_variants(
            geom_segmented, n_starting_points=4
        )
        for moved_start in starting_point_geoms:
            # Simplify strips noise below segment_length, but on straight edges
            # it also strips every densified vertex, leaving arbitrarily long
            # segments. Chaikin cuts corners at 1/4 of each segment's length,
            # so re-segmentize afterwards to cap the rounding at segment_length
            # scale.
            moved_start = moved_start.simplify(
                tolerance=segment_length,
                preserve_topology=True,
            ).segmentize(segment_length * _CHAIKIN_SEGMENT_FACTOR)
            moved_start = cast(Polygon, moved_start)
            # No reversed pass: on closed rings Chaikin is direction-invariant
            # (each edge {a, b} yields the same cut points 0.75a + 0.25b and
            # 0.25a + 0.75b either way), so a reversed variant duplicates the
            # forward one bit-for-bit and only inflates the union below.
            #
            # Cap pre-union iterations at 2: the merged result is simplified
            # at segment_length / 5 below, which erases detail finer than
            # iteration 3+ adds (cuts beyond iteration 2 move the boundary by
            # less than that tolerance), while each extra iteration doubles
            # the vertex count entering the expensive union. The final
            # post-union pass still runs the full smooth_iterations.
            smoothed = _chaikin_corner_cutting(
                geom=moved_start,
                num_iterations=min(smooth_iterations, 2),
            )
            geom_iterations.append(smoothed)

    else:
        moved_start = geom_segmented.simplify(
            tolerance=segment_length,
            preserve_topology=True,
        ).segmentize(segment_length * _CHAIKIN_SEGMENT_FACTOR)
        moved_start = cast(LineString, moved_start)
        smoothed = _chaikin_corner_cutting(
            geom=moved_start,
            num_iterations=smooth_iterations,
        )
        geom_iterations.append(smoothed)

    if isinstance(geom, Polygon):
        # A self-intersecting Chaikin variant repairs (via make_valid) into a
        # polygon plus zero-area line debris; keep only the area so the union
        # below stays polygonal.
        geom_iterations = [_polygonal_only(make_valid(g)) for g in geom_iterations]

        dissolved = _combine_variants(geom_iterations)

        dissolved_poly = dissolved.simplify(
            tolerance=segment_length / 5,
            preserve_topology=True,
        )

        # If the union is a MultiPolygon, take the largest geometry
        if isinstance(dissolved_poly, MultiPolygon):
            largest_geom = max(dissolved_poly.geoms, key=lambda x: x.area)
            dissolved_poly = largest_geom
    else:
        # LineStrings: skip make_valid/unary_union — self-intersecting lines
        # are geometrically valid, and unary_union would split them at
        # crossing points into a MultiLineString
        dissolved_poly = geom_iterations[0].simplify(
            tolerance=segment_length / 5,
            preserve_topology=True,
        )

    assert isinstance(dissolved_poly, (Polygon, LineString)), (
        f"Resulting geometry must be Polygon or LineString. Got {type(dissolved_poly)}."
    )

    # The simplify above can again leave long segments on straight stretches
    # (a Chaikin-rounded corner deviates < segment_length / 5 from the sharp
    # original, so simplify removes it entirely); re-segmentize so the final
    # Chaikin pass cannot re-round at a larger scale.
    dissolved_poly = dissolved_poly.segmentize(segment_length * _CHAIKIN_SEGMENT_FACTOR)

    smoothed_geom = _chaikin_corner_cutting(
        geom=dissolved_poly,
        num_iterations=smooth_iterations,
    )

    def finish(candidate: Polygon | LineString) -> Polygon | LineString:
        """Apply the optional area-preservation step to a smoothed result."""
        if (
            original_area is not None
            and isinstance(candidate, Polygon)
            and preserve_area
        ):
            # Convert percentage to absolute tolerance based on original area
            # (e.g. area_tolerance 0.01 = 0.01% error = 99.99% preservation)
            absolute_tolerance = original_area * (area_tolerance / 100.0)
            return _preserve_area_with_buffer(
                polygon=candidate,
                target_area=original_area,
                tolerance=absolute_tolerance,
            )
        return candidate

    smoothed_geom = finish(smoothed_geom)

    # Thin features fold in two ways that both surface as a sharp concave turn:
    # variants disagree about an arm ~one segment_length wide and their union
    # carries a forked slit, and — more severely — smoothing collapses such an
    # arm (it can shed ~40% of its area), so the area-preservation buffer above
    # has to expand it a long way and a large uniform outward buffer sharpens
    # the inner corner where arms meet into a fold/cusp. A correctly smoothed
    # result is gently curved on the concave side (thin arms may end in
    # legitimate sharp convex hairpins), so on a sharp concave turn fill that
    # notch with a small closing (dilate-erode), re-smooth, and restore area
    # again. Seal the area-preserved result rather than the pre-buffer union:
    # by now it already sits at target area, so this second area pass barely
    # moves the boundary and cannot re-introduce the fold the way sealing the
    # area-deficient union would once finish() buffered it back out.
    if (
        isinstance(smoothed_geom, Polygon)
        and _max_concave_turn_degrees(smoothed_geom) > 60
    ):
        radius = segment_length / 4
        sealed = make_valid(smoothed_geom.buffer(radius).buffer(-radius)).simplify(
            tolerance=segment_length / 5,
            preserve_topology=True,
        )
        if isinstance(sealed, MultiPolygon):
            sealed = max(sealed.geoms, key=lambda x: x.area)
        if isinstance(sealed, Polygon):
            resmoothed = _chaikin_corner_cutting(
                geom=sealed.segmentize(segment_length * _CHAIKIN_SEGMENT_FACTOR),
                num_iterations=smooth_iterations,
            )
            smoothed_geom = finish(resmoothed)

    return smoothed_geom
