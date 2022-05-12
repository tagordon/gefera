from setuptools import setup 

import pathlib
HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(name='gefera', 
      version='0.1',
      description='Gaussian processes for multi-wavelength and spectral observations',
      long_description=README,
      long_description_content_type="text/markdown",
      url='http://github.com/tagordon/specgp',
      author='Tyler Gordon',
      author_email='tagordon@uw.edu', 
      license='MIT',
      packages=['gefera'],
      install_requires=['numpy',
                        'scipy',
			'astropy'],
      zip_safe=False)
