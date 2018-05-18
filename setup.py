from setuptools import setup

setup(name='pycartogram',
      version='0.0.3',
      description='Collection of tools to effeciently generate cartograms as developed by Gastner and Newman.',
      url='https://github.com/benmaier/pycartogram',
      author='Benjamin F. Maier',
      author_email='benjaminfrankmaier@gmail.com',
      license='MIT',
      packages=['pycartogram'],
      include_package_data = True,
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'shapely',
          'descartes',
          'cartopy',
          'progressbar2',
          'visvalingamwyatt',
      ],
      dependency_links=[
          ],
      zip_safe=False)
