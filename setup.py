desc = '''
Elasto-plastic material model that usage a manifold of quadratic potentials.
'''

from setuptools import setup, Extension

import sys,re
import setuptools
import pybind11
import cppmat

header = open('src/ElastoPlasticQPot3d/ElastoPlasticQPot3d.h','r').read()
world  = re.split(r'(.*)(\#define ELASTOPLASTICQPOT3D_WORLD_VERSION\ )([0-9]+)(.*)',header)[3]
major  = re.split(r'(.*)(\#define ELASTOPLASTICQPOT3D_MAJOR_VERSION\ )([0-9]+)(.*)',header)[3]
minor  = re.split(r'(.*)(\#define ELASTOPLASTICQPOT3D_MINOR_VERSION\ )([0-9]+)(.*)',header)[3]

__version__ = '.'.join([world,major,minor])

ext_modules = [
  Extension(
    'ElastoPlasticQPot3d',
    ['src/ElastoPlasticQPot3d/python.cpp'],
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
      cppmat  .get_include(False),
      cppmat  .get_include(True ),
      cppmat  .find_eigen()
    ],
    language='c++'
  ),
]

setup(
  name             = 'ElastoPlasticQPot3d',
  description      = 'Elasto-plastic material model',
  long_description = desc,
  version          = __version__,
  license          = 'MIT',
  author           = 'Tom de Geus',
  author_email     = 'tom@geus.me',
  url              = 'https://github.com/tdegeus/ElastoPlasticQPot3d',
  ext_modules      = ext_modules,
  install_requires = ['pybind11>=2.2.0','cppmat>=0.4.1'],
  cmdclass         = {'build_ext': cppmat.BuildExt},
  zip_safe         = False,
)
