import numpy
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "qGI python module",
    ext_modules = cythonize('gimp.pyx'), # accepts a glob pattern
    include_dirs=[numpy.get_include()]
)
