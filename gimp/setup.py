import numpy
from distutils.core import setup
from Cython.Build import cythonize

setup(
      name = "qGI python module",
      version='0.1',
      description='',
      author='Matt Wittmann',
      author_email='mwittman@ucsc.edu',
      url='http://people.ucsc.edu/~mwittman',
      ext_modules = cythonize('gimp.pyx'), # accepts a glob pattern
      include_dirs=[numpy.get_include()]
     )
