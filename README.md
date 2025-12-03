# pycartogram

A Python package for generating cartograms using the diffusion method by
[Gastner and Newman (2004)](http://www.pnas.org/cgi/content/abstract/101/20/7499).

This package requires [cCartogram](https://github.com/benmaier/cCartogram), the Python port of Mark Newman's [original C code](http://www-personal.umich.edu/~mejn/cart/).

## Installation

```bash
# Install core package
pip install pycartogram

# Install with geographic extras (cartopy, geopandas)
pip install pycartogram[geo]

# Install development dependencies
pip install pycartogram[dev]
```

### Prerequisites

- **cCartogram**: Install from [GitHub](https://github.com/benmaier/cCartogram)
- **Cartopy** (optional, for `[geo]`): Requires system libraries `geos` and `proj`

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

2. Run `make html` in `docs/` — notebooks are automatically copied and converted to HTML by nbsphinx.

More examples can also be found in `sandbox/`.

## License

MIT
