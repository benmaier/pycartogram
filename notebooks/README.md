# Example Notebooks

Place Jupyter notebooks (`.ipynb` files) in this directory.

They will be automatically:
1. Copied to `docs/examples/` during documentation build
2. Converted to HTML by nbsphinx
3. Included in the "Examples" section of the documentation

## Naming Convention

Name notebooks with a numeric prefix for ordering:
- `01_basic_usage.ipynb`
- `02_point_cartograms.ipynb`
- `03_geopandas_integration.ipynb`

## Building Documentation

```bash
cd docs
make html
```

The Makefile will copy notebooks before building.
