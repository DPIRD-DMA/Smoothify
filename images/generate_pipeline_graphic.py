"""Generate the README pipeline graphic: how a polygon moves through smoothify.

Runs a small pixelated polygon through a simplified version of the real
smoothing pipeline (using the library's own building blocks, minus edge-case
handling) and renders every intermediate product as one multi-panel figure.

Usage:
    python images/generate_pipeline_graphic.py
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely import make_valid
from shapely.geometry import Polygon
from shapely.ops import unary_union

from smoothify.smoothify_core import (
    _CHAIKIN_SEGMENT_FACTOR,
    _chaikin_corner_cutting,
    _combine_variants,
    _generate_starting_point_variants,
    _polygonal_only,
    _preserve_area_with_buffer,
)

OUT = Path(__file__).parent / "pipeline_steps.png"
SEGMENT_LENGTH = 1.0  # "pixel size" of the demo shape
VARIANT_COLORS = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
# distinct linestyles so variants stay individually visible even where
# their outlines coincide
VARIANT_STYLES = ["-", "--", "-.", ":"]


def pixel_blob() -> Polygon:
    """A pond-like raster blob: a smooth shape rasterized at pixel size 1.

    Rasterizing a real curve gives an authentic staircase boundary at a
    realistic feature-to-pixel ratio (~20 pixels across). The lobe placement
    is chosen so all four start-point variants simplify to visibly different
    polygons (rotated anchors can otherwise land where two variants come out
    identical and overprint in the figure)."""
    xs, ys = np.meshgrid(np.arange(26), np.arange(20))
    cx, cy = xs + 0.5, ys + 0.5
    lobe_a = ((cx - 8.5) / 8.0) ** 2 + ((cy - 8) / 6.0) ** 2 < 1
    lobe_b = ((cx - 16) / 6.0) ** 2 + ((cy - 10.5) / 4.5) ** 2 < 1
    grid = lobe_a | lobe_b
    squares = [
        Polygon([(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)])
        for y, x in zip(*np.nonzero(grid), strict=True)
    ]
    merged = unary_union(squares)
    assert isinstance(merged, Polygon)
    # unary_union leaves a vertex at every pixel boundary along straight runs;
    # drop the collinear ones (tolerance 0) so only the true staircase corners
    # remain, matching what a polygonized raster's meaningful vertices are.
    corners = merged.simplify(0)
    assert isinstance(corners, Polygon)
    return corners


def draw(ax, geom, color="black", lw=1.8, dots=False, alpha=1.0, style="-"):
    xs, ys = geom.exterior.xy
    ax.plot(
        xs,
        ys,
        style,
        color=color,
        linewidth=lw,
        alpha=alpha,
        solid_capstyle="round",
    )
    if dots:
        ax.plot(xs, ys, "o", color=color, markersize=3.5, alpha=alpha)


def reference(ax, geom):
    gpd.GeoSeries([geom]).plot(ax=ax, color="0.92")
    draw(ax, geom, color="0.75", lw=1.0)


def main() -> None:
    original = pixel_blob()

    # --- the simplified pipeline, capturing intermediates -------------------
    # 1. variants that start at evenly spaced arc-length positions, each
    #    simplified (noise removal) and re-segmentized so corner rounding stays
    #    capped at segment_length. No up-front densify is needed: the start
    #    points are interpolated directly at their arc-length positions.
    variants = []
    for variant in _generate_starting_point_variants(original, n_starting_points=4):
        variant = variant.simplify(
            tolerance=SEGMENT_LENGTH, preserve_topology=True
        ).segmentize(SEGMENT_LENGTH * _CHAIKIN_SEGMENT_FACTOR)
        variants.append(variant)

    # 2. Chaikin corner cutting per variant
    smoothed_variants = [
        _polygonal_only(make_valid(_chaikin_corner_cutting(v, num_iterations=2)))
        for v in variants
    ]

    # 3. per-point median merge: a start-invariant consensus of the variants
    merged = _combine_variants(smoothed_variants).simplify(
        tolerance=SEGMENT_LENGTH / 5, preserve_topology=True
    )
    assert isinstance(merged, Polygon)

    # 4. final smoothing pass
    final_smooth = _chaikin_corner_cutting(
        merged.segmentize(SEGMENT_LENGTH * _CHAIKIN_SEGMENT_FACTOR), num_iterations=3
    )

    # 5. restore the original area by buffering
    final = _preserve_area_with_buffer(
        final_smooth, target_area=original.area, tolerance=original.area * 1e-4
    )

    # --- render --------------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2))
    panels = axes.ravel()

    ax = panels[0]
    gpd.GeoSeries([original]).plot(ax=ax, color="#cfe3f5")
    draw(ax, original, dots=True)
    ax.set_title("1. Pixelated input\n(polygonized raster)")

    ax = panels[1]
    reference(ax, original)
    for v, c, s in zip(variants, VARIANT_COLORS, VARIANT_STYLES, strict=True):
        draw(ax, v, color=c, lw=1.6, alpha=0.9, style=s)
    ax.set_title("2. Start at 4 arc-length\npositions, simplify each")

    ax = panels[2]
    reference(ax, original)
    for v, c, s in zip(smoothed_variants, VARIANT_COLORS, VARIANT_STYLES, strict=True):
        draw(ax, v, color=c, lw=1.6, alpha=0.9, style=s)
    ax.set_title("3. Chaikin corner cutting\nper variant")

    ax = panels[3]
    reference(ax, original)
    for v, c, s in zip(smoothed_variants, VARIANT_COLORS, VARIANT_STYLES, strict=True):
        draw(ax, v, color=c, lw=1.0, alpha=0.35, style=s)
    draw(ax, merged, color="black", lw=2.0)
    ax.set_title("4. Per-point median merge\n(start-invariant consensus)")

    ax = panels[4]
    reference(ax, original)
    draw(ax, final_smooth, color="black", lw=1.8)
    ax.set_title("5. Final smoothing pass")

    ax = panels[5]
    reference(ax, original)
    gpd.GeoSeries([final]).plot(ax=ax, color="#d3eed3", alpha=0.8)
    draw(ax, final, color="#1a7a1a", lw=2.0)
    err = abs(final.area - original.area) / original.area
    ax.set_title(f"6. Restore original area\n(area error {err:.4%})")

    minx, miny, maxx, maxy = original.buffer(1.4).bounds
    for ax in panels:
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(minx, maxx)
        ax.set_ylim(miny, maxy)
        for spine in ax.spines.values():
            spine.set_color("0.8")

    fig.suptitle("How smoothify works", fontsize=15, y=0.99)
    plt.tight_layout()
    plt.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
