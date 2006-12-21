from distutils.core import setup, Extension

_hfkmodule = Extension('hfk._hfk',
                       sources = ['_hfkmodule.cpp'],
                       include_dirs = ['.'])

setup (name = '_hfk',
       version = '1.0',
       description = 'Computes Heegaard Floer homology for links',
       packages = ['hfk'],
       ext_modules = [_hfkmodule])
