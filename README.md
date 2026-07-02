<p align="center">
  <img src="https://raw.githubusercontent.com/DPIRD-DMA/Smoothify/main/images/smoothify_logo.png" alt="Smoothify Text" width="600">
</p>

<p align="center">
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python Version"></a>
  <a href="https://github.com/DPIRD-DMA/Smoothify/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="License"></a>
  <a href="https://pypi.org/project/smoothify/"><img src="https://img.shields.io/pypi/v/smoothify.svg" alt="PyPI version"></a>
  <a href="https://pepy.tech/projects/smoothify"><img src="https://static.pepy.tech/badge/smoothify" alt="PyPI downloads"></a>
  <a href="https://anaconda.org/conda-forge/smoothify"><img src="https://img.shields.io/conda/vn/conda-forge/smoothify.svg" alt="Conda version"></a>
  <a href="https://github.com/DPIRD-DMA/Smoothify/tree/main/examples"><img src="https://img.shields.io/badge/Tutorials-Learn-brightgreen" alt="Tutorials"></a>
</p>

📋 [View Changelog](https://github.com/DPIRD-DMA/Smoothify/blob/main/CHANGELOG.md)

A Python package for smoothing and refining geometries derived from raster data classifications. Smoothify transforms jagged polygons and lines resulting from raster-to-vector conversion into smooth, visually appealing features using an optimized implementation of Chaikin's corner-cutting algorithm.

## Problem

Polygons and lines derived from classified raster data (e.g., ML model predictions, spectral indices, or remote sensing classifications) often have unnatural "stair-stepped" or "pixelated" edges that:

- Are visually unappealing in maps and GIS applications
- Can be difficult to work with in downstream vector processing
- Don't represent the real-world features they're meant to depict

## Solution

Smoothify applies an optimized implementation of Chaikin's corner-cutting algorithm along with other geometric processing to create smooth, natural-looking features while:

- Preserving the general shape and area of polygons
- Supporting all shapely geometry types
- Handling shapes with interior holes
- Efficiently processing large datasets with multiprocessing

<p align="left">
  <img src="https://raw.githubusercontent.com/DPIRD-DMA/Smoothify/main/images/smoothify_hero.png" alt="Smoothify Hero Image" width="600">
</p>

## Installation

```bash
uv add smoothify
```
or 
```bash
pip install smoothify
```

or 
```bash
conda install conda-forge::smoothify
```

## Quick Start

```python
import geopandas as gpd
from smoothify import smoothify

# Load your polygonized raster data
polygon_gdf = gpd.read_file("path/to/your/polygons.gpkg")

# Apply smoothing (segment_length auto-detected from geometry)
smoothed_gdf = smoothify(
    geom=polygon_gdf,
    smooth_iterations=3,  # More iterations = smoother result
    num_cores=4  # Use parallel processing for large datasets
)

# Or specify segment_length explicitly (generally recommended)
smoothed_gdf = smoothify(
    geom=polygon_gdf,
    segment_length=10.0,  # Use the original raster resolution
    smooth_iterations=3,
    num_cores=4
)

# Save the result
smoothed_gdf.to_file("smoothed_polygons.gpkg")
```

## Examples

Example notebooks:
- [Usage examples](https://github.com/DPIRD-DMA/Smoothify/blob/main/examples/usage_examples.ipynb)
- [Smoothify vs. Shapely comparison](https://github.com/DPIRD-DMA/Smoothify/blob/main/examples/smoothify_vs_shapely_comparison.ipynb)
- [Real-world water example](https://github.com/DPIRD-DMA/Smoothify/blob/main/examples/real_world_water_example.ipynb)
- [Merging holes](https://github.com/DPIRD-DMA/Smoothify/blob/main/examples/merge_holes_examples.ipynb)


### Basic Polygon Smoothing

Transform pixelated polygons from raster data into smooth, natural-looking features:

<p align="left">
  <img src="https://raw.githubusercontent.com/DPIRD-DMA/Smoothify/main/images/example_1_polygon.png" alt="Basic Polygon Smoothing" width="600">
</p>

### LineString Smoothing

Works perfectly for roads, streams, and other linear features:

<p align="left">
  <img src="https://raw.githubusercontent.com/DPIRD-DMA/Smoothify/main/images/example_2_linestring.png" alt="LineString Smoothing" width="600">
</p>

### Controlling Smoothness with Iterations

The `smooth_iterations` parameter controls how smooth the result will be:

<p align="left">
  <img src="https://raw.githubusercontent.com/DPIRD-DMA/Smoothify/main/images/example_3_iterations.png" alt="Effect of Different Iterations" width="700">
</p>

### Merging Adjacent Geometries

When processing multiple adjacent polygons, allowing merge_collection = True produces a combined result:

<p align="left">
  <img src="https://raw.githubusercontent.com/DPIRD-DMA/Smoothify/main/images/example_4_merging.png" alt="Merging Adjacent Geometries" width="700">
</p>


## General Usage

The `smoothify()` function accepts three types of input:

### 1. GeoDataFrame
```python
import geopandas as gpd
from smoothify import smoothify
# By default this will dissolve adjacent polygons before smoothing
gdf = gpd.read_file("polygons.gpkg")
smoothed_gdf = smoothify(
    geom=gdf,
    segment_length=10.0,
    smooth_iterations=3,
    num_cores=4
)

# Dissolve geometries by a specific field before smoothing
# Useful for merging adjacent polygons with the same classification
gdf_with_classes = gpd.read_file("classified_polygons.gpkg")
smoothed_by_class = smoothify(
    geom=gdf_with_classes,
    segment_length=10.0,
    smooth_iterations=3,
    merge_collection=True,
    merge_field="land_type",  # Merge adjacent geometries with same land_type
    num_cores=4
)
```

### 2. Single Geometry
```python
from shapely.geometry import Polygon
from smoothify import smoothify

polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
smoothed_polygon = smoothify(
    geom=polygon,
    smooth_iterations=3
)
```

### 3. List of Geometries or GeometryCollection
```python
from shapely.geometry import Polygon, LineString
from smoothify import smoothify

geometries = [
    Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
    LineString([(0, 0), (5, 5), (10, 0)])
]
smoothed = smoothify(
    geom=geometries,
    segment_length=1.0,
    smooth_iterations=3
)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `geom` | GeoDataFrame, BaseGeometry, or list[BaseGeometry] | Required | The geometry/geometries to smooth |
| `segment_length` | float | None | Resolution of the original raster data in map units. If None (default), automatically detects by finding the minimum segment length (from a data sample). Recommended to specify explicitly when known |
| `smooth_iterations` | int | 3 | Number of Chaikin corner-cutting iterations (typically 3-5). Higher values = smoother output with more vertices |
| `num_cores` | int | 0 | Number of CPU cores for parallel processing (0 = all available cores, 1 = serial) |
| `merge_collection` | bool | True | Whether to merge/dissolve adjacent geometries in collections before smoothing |
| `merge_field` | str | None | **GeoDataFrame only**: Column name to use for dissolving geometries. Only valid when `merge_collection=True`. If None, dissolves all geometries together. If specified, dissolves geometries grouped by the column values |
| `merge_multipolygons` | bool | True | Whether to merge adjacent polygons within MultiPolygons before smoothing |
| `merge_holes` | bool | True | Whether to join holes that touch or nearly touch (e.g. diagonally adjacent raster cells) before smoothing, so they smooth into one coherent opening instead of separate overlapping shapes |
| `preserve_area` | bool | True | Whether to restore original area after smoothing via buffering (applies to Polygons only) |
| `area_tolerance` | float | 0.01 | Percentage of original area allowed as error (e.g., 0.01 = 0.01% error = 99.99% preservation). Only affects Polygons when preserve_area=True |

## How It Works

Smoothify uses an advanced multi-step smoothing pipeline. The numbered steps below correspond to the panels in the figure:

<p align="left">
  <img src="https://raw.githubusercontent.com/DPIRD-DMA/Smoothify/main/images/pipeline_steps.png" alt="Smoothify pipeline steps" width="800">
</p>


1. **Pixelated input** — a polygon straight from raster-to-vector conversion, with a stair-stepped boundary
2. **Multiple variants** (for Polygons) that start at evenly spaced arc-length positions, each simplified to strip staircase noise, so no artifact is tied to a fixed start vertex
3. **Chaikin corner cutting** applied to each variant
4. **Per-point median merge** — a start-invariant consensus that resolves the variants' disagreements
5. **Final smoothing pass** on the merged result
6. **Restore original area** via buffering (for Polygons, when `preserve_area=True`)

Two steps are not shown in the figure: before step 2, touching holes are joined (for Polygons, when `merge_holes=True`) so they smooth as one opening; and after step 6, the result is checked for sharp concave folds left by features near the smoothing scale (e.g. one-pixel-wide arms) and repaired with a small morphological opening/closing bounded at `segment_length / 4`.

## Invalid Geometries

Smoothify does not repair invalid input. If it encounters an invalid geometry (e.g. a self-intersecting polygon), it returns that geometry **unchanged** and emits a warning, instead of crashing or silently producing an empty geometry. This is consistent whether you pass a single geometry, a list/collection, or a GeoDataFrame.

If you want invalid geometries smoothed, repair them first with shapely's `make_valid()`:

```python
# GeoDataFrame
gdf.geometry = gdf.geometry.make_valid()
smoothed_gdf = smoothify(gdf, segment_length=10.0)

# Single geometry
from shapely import make_valid
smoothed = smoothify(make_valid(polygon), segment_length=1.0)
```

## Performance Considerations

- **Parallel Processing**: For large GeoDataFrames or collections, use `num_cores` = 0 to enable parallel processing
- **Duplicate Shapes**: Geometries that are translated copies of the same shape (common in raster-derived data, e.g. single-pixel polygons) are automatically smoothed once and the result reused
- **Smoothing Iterations**: Values of 3-5 typically provide good results. Higher values create smoother output but increase processing time and vertex count
- **Memory Usage**: Scales with geometry complexity. The algorithm creates multiple variants during smoothing
- **Optimal segment_length**: Anything from about half the original raster pixel size and up should produce reasonable output — larger values produce more rounded output, smaller values stay more faithful to the original geometry

## Running the Tests

Smoothify uses [pytest](https://pytest.org/). After cloning the repository, install the development dependencies and run the suite with [uv](https://docs.astral.sh/uv/):

```bash
# Install dependencies (including the dev group)
uv sync

# Run all tests (parallelised across CPU cores by default via pytest-xdist)
uv run pytest tests/

# Run with coverage
uv run pytest tests/ --cov=smoothify --cov-report=html

# Run a single test (add `-n 0` to disable parallelism for clearer output)
uv run pytest tests/test_chaikin.py::TestChaikinCornerCutting::test_simple_square_polygon -n 0
```

The suite runs in parallel by default (`-n auto` in `pytest.ini`); pass `-n 0` to run serially when debugging. If you prefer not to use uv, install the dev dependencies into your environment and run `pytest tests/` directly.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/DPIRD-DMA/Smoothify/blob/main/LICENSE) file for details.
