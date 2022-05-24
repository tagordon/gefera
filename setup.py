# from setuptools import setup, Extension
from numpy.distutils.core import setup, Extension

import pathlib
HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

keplib = Extension(name = 'gefera.keplib',
                   sources = ['fortran/kep.f90']
                  )

hrchlib = Extension(name = 'gefera.hrchlib', 
                    sources = ['fortran/kep.f90', 
                               'fortran/hrch.f90'
                              ]
                   )
conflib = Extension(name = 'gefera.conflib', 
                    sources = ['fortran/kep.f90', 
                               'fortran/conf.f90'
                              ]
                   )
photlib = Extension(name = 'gefera.photlib', 
                    sources = ['fortran/ellip.f90', 
                               'fortran/phot_nograd.f90',
                               'fortran/phot.f90'
                              ]
                   )

setup(name='gefera', 
      version='0.1',
      description='Light curves for mutual transits of limb-darkened stars',
      long_description=README,
      long_description_content_type="text/markdown",
      url='http://github.com/tagordon/gefera',
      author='Tyler Gordon',
      author_email='tagordon@uw.edu', 
      license='MIT',
      packages=['gefera'],
      install_requires=['numpy',
                        'scipy',
                        'astropy'
                       ],
      zip_safe=True,
      ext_modules = [hrchlib, 
                     conflib, 
                     photlib,
                     keplib
                    ]
      )
