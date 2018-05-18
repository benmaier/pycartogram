# pycartogram

A Python package providing a bunch of tools to procude cartograms based on the method by
[Gastner and Newman (2004)](http://www.pnas.org/cgi/content/abstract/101/20/7499]).

This package does not work without [cCartogram](https://github.com/benmaier/cCartogram) which is the Python port of Mark Newman's [original code](http://www-personal.umich.edu/~mejn/cart/).


## Examples

Find an example below and more in subfolder `sandbox`.

```python
import matplotlib.pyplot as pl
from shapely.geometry import Polygon
from pycartogram import WardCartogram

A = Polygon([ (0.,0.), (1.,0.), (1.,1.), (0.,1.), (0.,0.), ])
B = Polygon([ (1.,0.), (2.,0.), (2.,1.), (1.,1.), (1.,0.), ])
C = Polygon([ (2.,0.), (3.,0.), (3.,1.), (2.,1.), (2.,0.), ])

wards = [A,B,C]
ward_population = [2.,20.,200.]

# create cartogram instance from wards
carto = WardCartogram(wards=wards,
                      ward_density=ward_population,
                      norm_density=True, # needs to be normed by area
                      margin_ratio=0.5,
                      x_raster_size=128,
                      y_raster_size=64,)

# compute density matrix and cartogram, and transform
# polygon coordinates accordingly
carto.compute(verbose=True)

# plot old wards color
carto.plot(show_new_wards=False,
           ward_colors = [0.3,0.5,0.8],
           edge_colors = 'k',
           bg_color = 'k',
          )

# plot new wards with old density
fig, ax = carto.plot(show_new_wards=True,
           edge_colors = 'k'
          )

pl.show()
```

## Install

For all systems, first clone this repository.

Be sure you have [cCartogram](https://github.com/benmaier/cCartogram) installed.
Even though `pip` should take care of it, `cCartogram` needs some extra libraries such that 
`pip` might fail.

Another package needed is [Cartopy](https://github.com/SciTools/cartopy) which requires the C-libraries `geos` and `proj`, which I installed using MacPorts. Installing Cartopy with MacPorts failed during building, which is, I think, related to the fact that the library and include directories of `proj` are not found automatically. I ended up cloning the Cartopy source code repo and manually including the paths `/opt/local/lib/proj49/lib` and `/opt/local/lib/proj49/include` for installation.

To install do

    $ sudo pip install ./pycartogram
