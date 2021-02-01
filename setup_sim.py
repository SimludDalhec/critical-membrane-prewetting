#!/usr/bin/python
import numpy
import pyximport; pyximport.install()
from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import sys

setup(name = 'SetupTb',
      ext_modules=cythonize(['ising_class.pyx','bulk_polymer.pyx','tether.pyx','polymer_membrane_sim.pyx']),
      include_dir=[numpy.get_include()],
      cmdclass={'build_ext': build_ext}
)

