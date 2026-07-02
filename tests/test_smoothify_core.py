"""Tests for core smoothify geometry functions."""

import numpy as np
import pytest
from shapely.geometry import LineString, Polygon
from shapely.geometry.base import BaseGeometry

from smoothify import smoothify_core
from smoothify.smoothify_core import (
    _generate_starting_point_variants,
    _join_adjacent,
    _preserve_area_brentq,
    _preserve_area_with_buffer,
    _roll0,
    _rotate_polygon_start,
    _rotate_ring_coords,
    _smoothify_geometry,
)


class TestRotatePolygonStart:
    """Test suite for polygon rotation functions."""

    def test_rotate_square(self):
        """Test rotating a square's starting point."""
        square = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        rotated = _rotate_polygon_start(square, shift=0.25)

        assert isinstance(rotated, Polygon)
        assert rotated.is_valid
        # Area should be identical
        assert abs(rotated.area - square.area) < 1e-10

    def test_rotate_preserves_geometry(self):
        """Rotation must yield the identical shape, only a different start.

        The arc-length start vertex is interpolated on an existing edge, so it
        is collinear and the ring is geometrically unchanged. Uses an irregular
        polygon (asymmetric, uneven vertex spacing) so a placement bug would
        actually move the boundary rather than hide behind symmetry."""
        poly = Polygon([(0, 0), (30, 0), (30, 10), (18, 10), (18, 4), (0, 4)])
        for shift in (0.1, 0.25, 0.5, 0.73, 0.9):
            rotated = _rotate_polygon_start(poly, shift=shift)
            assert rotated.is_valid
            assert abs(rotated.area - poly.area) < 1e-9
            # Same footprint: symmetric difference is (near) zero area.
            assert poly.symmetric_difference(rotated).area < 1e-9

    def test_rotate_start_at_arc_length(self):
        """The new start sits at the requested arc-length fraction, not a vertex
        index. On this 20x10 rectangle (perimeter 60) shift=0.25 is 15 units
        along the first (length-20) edge, i.e. the interpolated point (15, 0) -
        whereas index-based rotation would have jumped to the corner (20, 0)."""
        rect = Polygon([(0, 0), (20, 0), (20, 10), (0, 10)])
        coords = _rotate_ring_coords(rect.exterior, shift=0.25)
        assert coords[0] == pytest.approx((15.0, 0.0))
        # Ring stays closed.
        assert coords[0] == pytest.approx(coords[-1])


class TestGenerateStartingPointVariants:
    """Test suite for generating polygon variants."""

    def test_polygon_variants(self):
        """Test generating variants for a polygon."""
        polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        variants = _generate_starting_point_variants(polygon, n_starting_points=4)

        assert len(variants) == 4
        assert all(isinstance(v, Polygon) for v in variants)
        assert all(v.is_valid for v in variants)

    def test_linestring_variants(self):
        """Test that linestrings return single variant."""
        line = LineString([(0, 0), (5, 5), (10, 0)])
        variants = _generate_starting_point_variants(line, n_starting_points=4)

        # LineStrings should return single item
        assert len(variants) == 1
        assert variants[0] == line


class TestPreserveAreaWithBuffer:
    """Test suite for area preservation function."""

    def test_preserve_area_exact(self):
        """Test preserving area when already exact."""
        polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        original_area = polygon.area

        preserved = _preserve_area_with_buffer(
            polygon, target_area=original_area, tolerance=1e-6
        )

        # Should return same polygon if area already matches
        assert abs(preserved.area - original_area) < 1e-6

    def test_preserve_area_smaller_polygon(self):
        """Test expanding a smaller polygon to target area."""
        small_polygon = Polygon([(0, 0), (5, 0), (5, 5), (0, 5)])
        target_area = 100.0  # Larger than 25

        preserved = _preserve_area_with_buffer(
            small_polygon, target_area=target_area, tolerance=1e-3
        )

        # Should be close to target area
        assert abs(preserved.area - target_area) < 1e-3
        assert preserved.area > small_polygon.area

    def test_preserve_area_larger_polygon(self):
        """Test shrinking a larger polygon to target area."""
        large_polygon = Polygon([(0, 0), (20, 0), (20, 20), (0, 20)])
        target_area = 100.0  # Smaller than 400

        preserved = _preserve_area_with_buffer(
            large_polygon, target_area=target_area, tolerance=1e-3
        )

        # Should be close to target area
        assert abs(preserved.area - target_area) < 1e-3
        assert preserved.area < large_polygon.area

    def test_preserve_area_empty(self):
        """Empty input is returned unchanged."""
        empty = Polygon()
        assert _preserve_area_with_buffer(empty, target_area=10.0).is_empty

    def test_preserve_area_concave_shape(self):
        """Newton path reaches tolerance on a concave (non-convex) polygon."""
        # L-shape: the pi*d^2 Steiner seed assumes total turning of 2*pi, which
        # a reflex corner violates, so this exercises the Newton correction.
        l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
        for target in (l_shape.area * 1.05, l_shape.area * 0.95):
            preserved = _preserve_area_with_buffer(
                l_shape, target_area=target, tolerance=1e-4
            )
            assert abs(preserved.area - target) < 1e-4

    def test_preserve_area_uses_few_buffers(self):
        """The Newton solve should reach tolerance in far fewer buffers than the
        bracketed fallback would (guards the optimisation against regressions)."""
        polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        calls = {"n": 0}
        original_buffer = BaseGeometry.buffer

        def counting_buffer(self, *args, **kwargs):
            calls["n"] += 1
            return original_buffer(self, *args, **kwargs)

        BaseGeometry.buffer = counting_buffer
        try:
            preserved = _preserve_area_with_buffer(
                polygon, target_area=polygon.area * 1.1, tolerance=1e-4
            )
        finally:
            BaseGeometry.buffer = original_buffer

        assert abs(preserved.area - polygon.area * 1.1) < 1e-4
        assert calls["n"] <= 4

    def test_preserve_area_falls_back_to_brentq(self, monkeypatch):
        """When Newton is denied any steps the function must still reach
        tolerance via the bracketed Brent's-method fallback."""
        monkeypatch.setattr(smoothify_core, "_AREA_NEWTON_MAX_STEPS", 0)
        polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        for target in (polygon.area * 1.1, polygon.area * 0.9):
            preserved = _preserve_area_with_buffer(
                polygon, target_area=target, tolerance=1e-3
            )
            assert abs(preserved.area - target) < 1e-3

    def test_newton_and_brentq_agree(self):
        """The Newton path and the bracketed fallback land on the same area."""
        polygon = Polygon([(0, 0), (8, 0), (8, 8), (0, 8)])
        target = polygon.area * 1.07
        newton = _preserve_area_with_buffer(polygon, target_area=target, tolerance=1e-4)
        brentq = _preserve_area_brentq(
            polygon,
            target_area=target,
            tolerance=1e-4,
            current_area=polygon.area,
            perimeter=polygon.length,
        )
        assert brentq is not None
        assert abs(newton.area - target) < 1e-4
        assert abs(brentq.area - target) < 1e-4
        assert abs(newton.area - brentq.area) < 1e-3


class TestPreserveAreaBrentq:
    """Test suite for the bracketed Brent's-method fallback."""

    def test_brentq_grows_polygon(self):
        """Fallback expands a polygon to a larger target area."""
        polygon = Polygon([(0, 0), (5, 0), (5, 5), (0, 5)])
        result = _preserve_area_brentq(
            polygon,
            target_area=100.0,
            tolerance=1e-3,
            current_area=polygon.area,
            perimeter=polygon.length,
        )
        assert result is not None
        assert abs(result.area - 100.0) < 1e-3
        assert result.area > polygon.area

    def test_brentq_shrinks_polygon(self):
        """Fallback shrinks a polygon to a smaller target area."""
        polygon = Polygon([(0, 0), (20, 0), (20, 20), (0, 20)])
        result = _preserve_area_brentq(
            polygon,
            target_area=100.0,
            tolerance=1e-3,
            current_area=polygon.area,
            perimeter=polygon.length,
        )
        assert result is not None
        assert abs(result.area - 100.0) < 1e-3
        assert result.area < polygon.area


class TestJoinAdjacent:
    """Test suite for joining adjacent geometries."""

    def test_join_two_squares(self):
        """Test joining two adjacent squares."""
        square1 = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        square2 = Polygon([(10, 0), (20, 0), (20, 10), (10, 10)])

        merged = _join_adjacent([square1, square2], segment_length=10.0)

        assert merged.is_valid
        # Should create a single merged geometry
        assert isinstance(merged, Polygon) or hasattr(merged, "geoms")

    def test_join_single_geometry(self):
        """Test joining a single geometry."""
        polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        merged = _join_adjacent(polygon, segment_length=10.0)

        assert merged.is_valid


class TestSmoothifyGeometry:
    """Test suite for the main smoothify_geometry function."""

    def test_smoothify_simple_polygon(self):
        """Test smoothing a simple polygon."""
        square = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        smoothed = _smoothify_geometry(
            square, segment_length=1.0, smooth_iterations=3, preserve_area=True
        )

        assert isinstance(smoothed, Polygon)
        assert smoothed.is_valid
        # Should have more vertices
        assert len(smoothed.exterior.coords) > len(square.exterior.coords)

    def test_smoothify_linestring(self):
        """Test smoothing a linestring."""
        line = LineString([(0, 0), (5, 5), (10, 0)])
        smoothed = _smoothify_geometry(
            line, segment_length=1.0, smooth_iterations=3, preserve_area=False
        )

        assert isinstance(smoothed, LineString)
        assert smoothed.is_valid
        # Endpoints should be preserved
        assert smoothed.coords[0] == line.coords[0]
        assert smoothed.coords[-1] == line.coords[-1]

    def test_smoothify_self_intersecting_linestring(self):
        """Test that a self-intersecting LineString is smoothed without error.

        Lines that cross themselves are geometrically valid. The smoothing
        pipeline must not split them into a MultiLineString via make_valid
        or unary_union, which node lines at self-intersection points.
        """
        # S-curve that crosses itself — triggers unary_union splitting
        line = LineString(
            [
                (0, 0),
                (2, 0),
                (3, 0),
                (4, 1),
                (3, 2),
                (2, 2),
                (1, 1),
                (2, 0.5),
                (3, 0.5),
                (4, 0),
                (5, 0),
                (7, 0),
            ]
        )
        assert not line.is_simple  # confirm it self-intersects

        smoothed = _smoothify_geometry(
            line, segment_length=0.5, smooth_iterations=3, preserve_area=False
        )

        assert isinstance(smoothed, LineString)

    def test_preserve_area_option(self):
        """Test that preserve_area option works."""
        square = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
        original_area = square.area

        smoothed = _smoothify_geometry(
            square, segment_length=10.0, smooth_iterations=3, preserve_area=True
        )

        # Area should be close to original
        assert abs(smoothed.area - original_area) / original_area < 0.05  # Within 5%

    def test_no_preserve_area(self):
        """Test smoothing without area preservation."""
        square = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
        original_area = square.area

        smoothed = _smoothify_geometry(
            square, segment_length=10.0, smooth_iterations=3, preserve_area=False
        )

        # Area will likely be different (smaller due to corner cutting)
        assert smoothed.is_valid
        assert smoothed.area < original_area

    def test_invalid_geometry_type(self):
        """Test that invalid geometry type raises error."""
        with pytest.raises((ValueError, AttributeError)):
            _smoothify_geometry(
                "not a geometry", segment_length=1.0, smooth_iterations=3
            )  # type: ignore


class TestRoll0:
    """_roll0 is a hot-path stand-in for np.roll(a, k, axis=0); it must match
    it bit-for-bit for the shifts the pipeline uses (Chaikin's -1, and the
    +/-k alignment shifts)."""

    def test_matches_numpy_roll(self):
        rng = np.random.default_rng(0)
        for n in (1, 2, 5, 17, 240):
            a = rng.standard_normal((n, 2))
            for k in (-1, 0, 1, 3, -4, n, n + 2, -n - 1):
                np.testing.assert_array_equal(_roll0(a, k), np.roll(a, k, axis=0))

    def test_zero_shift_returns_same_object(self):
        """k == 0 (mod n) should short-circuit without copying."""
        a = np.arange(10, dtype=np.float64).reshape(5, 2)
        assert _roll0(a, 0) is a
        assert _roll0(a, 5) is a
