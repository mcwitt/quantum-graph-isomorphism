import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# from Cython.Build import cythonize

ext_modules = [Extension(
    name="qaa",
    sources=["qaa.pyx"],
        extra_objects=["qaa.o"],  # if you compile fc.cpp separately
    include_dirs = [
        '../include',
        numpy.get_include()
        ],  # .../site-packages/numpy/core/include
    language="c",
        # libraries=
        # extra_compile_args = "...".split(),
        # extra_link_args = "...".split()
    )]

setup(
    name = 'qaa',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
        # ext_modules = cythonize(ext_modules)  ? not in 0.14.1
    # version=
    # description=
    # author=
    # author_email=
    )

