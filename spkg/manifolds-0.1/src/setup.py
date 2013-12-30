from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import Cython.Compiler.Options
import os

Cython.Compiler.Options.old_style_globals = True

setup(name='manifolds',
	version='0.1',
	description='SageManifolds',
	author='Eric Gourgoulhon, Michal Bejger',
	author_email='eric.gourgoulhon@obspm.fr, bejger@camk.edu.pl',
	url='http://sagemanifolds.obspm.fr',
    license = "GPL v3",
	packages=['manifolds'],
    ext_modules=[
    Extension('manifolds.chart',
    sources = [os.path.join('manifolds','chart.pyx')]),
    ], 
    cmdclass = {'build_ext': build_ext}
)

