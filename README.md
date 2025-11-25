<p align="left">
  <img src="https://raw.githubusercontent.com/DPIRD-DMA/Smoothify/main/images/smoothify_logo.png" alt="Smoothify Text" width="600">
</p>

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/DPIRD-DMA/Smoothify/blob/main/LICENSE)
[![PyPI version](https://img.shields.io/pypi/v/smoothify.svg)](https://pypi.org/project/smoothify/)
[![Tutorials](https://img.shields.io/badge/Tutorials-Learn-brightgreen)](https://github.com/DPIRD-DMA/Smoothify/tree/main/examples)

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
- Supporting all shapley geometry types
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

## General Usage

The `smoothify()` function accepts three types of input:

### 1. GeoDataFrame
```python
import geopandas as gpd
from smoothify import smoothify

gdf = gpd.read_file("polygons.gpkg")
smoothed_gdf = smoothify(
    geom=gdf,
    segment_length=10.0,
    smooth_iterations=3,
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
| `merge_multipolygons` | bool | True | Whether to merge adjacent polygons within MultiPolygons before smoothing |
| `preserve_area` | bool | True | Whether to restore original area after smoothing via buffering (applies to Polygons only) |
| `area_tolerance` | float | 0.01 | Percentage of original area allowed as error (e.g., 0.01 = 0.01% error = 99.99% preservation). Only affects Polygons when preserve_area=True |

## How It Works

Smoothify uses an advanced multi-step smoothing pipeline:


1. Adds intermediate vertices along line segments (segmentize)
2. Generates multiple rotated variants (for Polygons) to avoid artifacts
3. Simplifies each variant to remove noise
4. Applies Chaikin corner cutting to smooth
5. Merges all variants via union to eliminate start-point artifacts
6. Applies final smoothing pass
7. Optionally restores original area via buffering (for Polygons)

## Performance Considerations

- **Parallel Processing**: For large GeoDataFrames or collections, use `num_cores` = 0 to enable parallel processing
- **Smoothing Iterations**: Values of 3-5 typically provide good results. Higher values create smoother output but increase processing time and vertex count
- **Memory Usage**: Scales with geometry complexity. The algorithm creates multiple variants during smoothing
- **Optimal segment_length**: Should match the original raster cell size (pixel size) or be slightly larger for best results

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/DPIRD-DMA/Smoothify/blob/main/LICENSE) file for details.