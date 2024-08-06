#!/usr/bin/env python

import sys
from distutils.core import setup, Extension

setup(name = 'MODYSA: MOlecular DYnamics Simulation Analysis',
      version = '0.1',
      description = 'A set of tools for analysis of MD simulation trajectories',
      author = 'Krzysztof Murzyn, Marcin Kurdziel',
      author_email = 'krzymu@gmail.com',
      packages = ['modysa', 'modysa.analyse', 'modysa.core', 'modysa.io', 
                  'modysa.io.gromacs', 'modysa.io.amber', 'modysa.misc',
                  'modysa.statistics', 'modysa.demo'],
      package_data = {'modysa': ['demo/data/*']}, 
      ext_package = 'modysa.'+sys.platform,
      ext_modules = [ Extension('gridwrapper', ['src/neighbors.c', 
	                                        'src/gridwrapper.c'], 
                                extra_compile_args=['-fno-strict-aliasing', 
                                                    '-DDOUBLE']),
                      Extension('caux', ['src/caux.c'])
                    ]
     )
