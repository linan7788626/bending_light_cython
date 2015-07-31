from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np
#import sys
#print sys.platform
#if sys.platform == "linux2" :
#    include_gsl_dir = "/usr/local/include/"
#    lib_gsl_dir = "/usr/local/lib/"
#elif sys.platform == "win32":
#    include_gsl_dir = r"c:\msys\1.0\local\include"
#    lib_gsl_dir = r"c:\msys\1.0\local\lib"

#ext = Extension("libv2_cv", ["./v2_cv.pyx","./cv_test.c","./fcic.c"],
#ext = Extension("libv3_cv", ["./v4_cv.pyx","./all_cv_test.c"],
ext = Extension("libfft_lensing", ["./fft_lensing.pyx","./lensing_funcs.c","./mycosmology.c"],
    include_dirs=[np.get_include(),
                  "./"],
    library_dirs=["./"],
    libraries=["m","fftw3","gsl","gslcblas"]
)

setup(ext_modules=[ext],
    cmdclass = {'build_ext': build_ext})
