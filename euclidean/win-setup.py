#Useful commands:

#Build a .pyd file
#python win-setup.py build --compiler=mingw32

#Build an installer
#python win-setup.py bdist_wininst

from distutils.core import setup, Extension
setup (name = "euclidean",
    version = "1.0",
    maintainer = "Michael Seiler",
    maintainer_email = "miseiler@gmail.com",
    description = "Euclidean C Ext for Windows",
    ext_modules = [Extension('euclidean',
                   extra_compile_args=['-O3', '-ffast-math'],
                   include_dirs=['/usr/lib/python2.5/site-packages/numpy/core/include/numpy', 'C:\\Python25\\Lib\\site-packages\\numpy-1.0.4.0002-py2.5-win32.egg\\numpy\\core\\include\\numpy'],
                   sources=['euclidean.c'])])
