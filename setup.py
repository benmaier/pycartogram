from setuptools import setup

setup(name='pycartogram',
      version='0.0.1',
      description='',
      url='https://github.com/benmaier/pycartogram',
      author='Benjamin F. Maier',
      author_email='benjaminfrankmaier@gmail.com',
      license='MIT',
      packages=['qsuite'],
      include_package_data = True,
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'shapely',
          'descartes',
          'cartopy',
          'progressbar',
      ],
      dependency_links=[
          ],
      zip_safe=False)
