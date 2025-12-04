# pycartogram

A Python package for generating cartograms using the diffusion method by
[Gastner and Newman (2004)](http://www.pnas.org/cgi/content/abstract/101/20/7499).

This package requires [cCartogram](https://github.com/benmaier/cCartogram), the Python port of Mark Newman's [original C code](http://www-personal.umich.edu/~mejn/cart/).

## Installation

Neither cCartogram nor pycartogram are on PyPI or conda-forge yet. Install from source in this order:

### 1. Install FFTW (required by cCartogram)

**macOS (Homebrew):**
```bash
brew install fftw
```

**Ubuntu/Debian:**
```bash
sudo apt-get install libfftw3-dev
```

**conda:**
```bash
conda install -c conda-forge fftw
```

### 2. Install cCartogram

```bash
pip install git+https://github.com/benmaier/cCartogram.git
```

### 3. Install pycartogram

```bash
# Core package
pip install git+https://github.com/benmaier/pycartogram.git

# With geographic extras (cartopy, geopandas)
pip install "pycartogram[geo] @ git+https://github.com/benmaier/pycartogram.git"
```

### Optional: Cartopy system dependencies

For the `[geo]` extras, Cartopy requires system libraries:

**macOS:**
```bash
brew install geos proj
```

**Ubuntu/Debian:**
```bash
sudo apt-get install libgeos-dev libproj-dev
```

### Troubleshooting: macOS LC_RPATH Error

If you encounter an error like:

```
ImportError: dlopen(.../cCartogram.cpython-312-darwin.so...): tried: '...' (duplicate LC_RPATH '/some/path/lib')
```

This happens when cCartogram's compiled library contains duplicate rpath entries. The error message shows the duplicate path. Fix it using `install_name_tool`:

```bash
# The .so file path is shown in the error message. Use it directly:
install_name_tool -delete_rpath '/the/duplicate/path/from/error' '/path/to/cCartogram.cpython-312-darwin.so'

# Example with a real path from the error:
install_name_tool -delete_rpath '/Users/you/miniconda3/lib' '/Users/you/miniconda3/lib/python3.12/site-packages/cCartogram.cpython-312-darwin.so'
```

You can verify the fix worked by checking the rpaths:

```bash
otool -l /path/to/cCartogram.cpython-312-darwin.so | grep -A2 LC_RPATH
```

Alternatively, reinstall cCartogram from source:

```bash
pip install --force-reinstall git+https://github.com/benmaier/cCartogram.git
```

## Features

- **WardCartogram**: Create cartograms from predefined ward/region boundaries with density values
- **PointCartogram**: Generate cartograms directly from point location data
- **VoronoiCartogram**: Create Voronoi tessellation-based cartograms
- **GeoDataFrameWardCartogram**: Direct GeoPandas integration with animated interpolation support
- **GoogleShapeProject**: Load and project Google location history data

### Utilities

- `fix_invalid_geometry()` / `fix_geodataframe_geometries()`: Repair self-intersecting polygons
- Smooth easing-based animation interpolation between original and cartogram geometries
- Export to JSON for web visualization

## Quick Start

```python
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from pycartogram import WardCartogram

# Define three adjacent rectangular wards
A = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
B = Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])
C = Polygon([(2, 0), (3, 0), (3, 1), (2, 1)])

wards = [A, B, C]
ward_population = [2., 20., 200.]

# Create and compute cartogram
carto = WardCartogram(
    wards=wards,
    ward_density=ward_population,
    norm_density=True,
    margin_ratio=0.5,
    x_raster_size=128,
    y_raster_size=64,
)
carto.compute(verbose=True)

# Plot result
fig, ax = carto.plot(show_new_wards=True, edge_colors='k')
plt.show()
```

![Quick Start Example](https://raw.githubusercontent.com/benmaier/pycartogram/main/sandbox/img/example_quickstart.png)

### With GeoPandas

```python
import geopandas as gpd
from pycartogram import GeoDataFrameWardCartogram

# Load your geodata
gdf = gpd.read_file("regions.geojson")
gdf['density'] = gdf['population'] / gdf.geometry.area

# Create cartogram
carto = GeoDataFrameWardCartogram(gdf, 'density')
carto.compute(verbose=True)

# Get transformed GeoDataFrame
new_gdf = carto.get_cartogram_geo_df()

# Or get interpolated frame for animation (t=0 original, t=1 cartogram)
interp_gdf = carto.get_interpolated_geo_df(t=0.5)
```

## Documentation

Build the documentation locally:

```bash
cd docs
pip install -r requirements.txt
make html
open _build/html/index.html
```

### Adding Examples

Example notebooks are stored in `notebooks/` and automatically included in the documentation:

1. Create a Jupyter notebook in `notebooks/` with a numeric prefix for ordering:
   - `01_basic_usage.ipynb`
   - `02_point_cartograms.ipynb`
   - `03_geopandas_integration.ipynb`

2. For notebooks with video output, place video files in `notebooks/animation/`

3. Run `make html` in `docs/` — notebooks and animations are automatically copied and converted to HTML by nbsphinx.

More examples can also be found in `sandbox/`.

## License

MIT
