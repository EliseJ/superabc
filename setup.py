#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
from pkg_resources import parse_version
from setuptools import setup, find_packages


numpy_min_version = '1.8'

def get_numpy_status():
	"""
	Returns a dictionary containing a boolean specifying whether NumPy
	is up-to-date, along with the version string (empty string if
	not installed).
	"""
	numpy_status = {}
	try:
		import numpy
		numpy_version = numpy.__version__
		numpy_status['up_to_date'] = parse_version(numpy_version) >= parse_version(numpy_min_version)
		numpy_status['version'] = numpy_version
	except ImportError:
		numpy_status['up_to_date'] = False
		numpy_status['version'] = ""
	return numpy_status

def setup_superabc():

	numpy_status = get_numpy_status()
	numpy_req_str = "superABC requires NumPy >= {0}.\n".format(numpy_min_version)
        
	if numpy_status['up_to_date'] is False:
		if numpy_status['version']:
            		raise ImportError("Your installation of NumPy""{0} is out-of-date.\n{1}".format(numpy_status['version'],numpy_req_str))
		else:
			raise ImportError("NumPy is not installed.\n{0}".format(numpy_req_str))
                              
    
	from numpy.distutils.misc_util import Configuration
	from numpy.distutils.core import setup
	
	setup(	name='superabc',
		version='1.1.0',
		author="Elise Jennings",
		author_email="elise.jennings@gmail.com ",
		url="https://github.com/EliseJ/superabc",
		description='A Python sampling method for obtaining cosmological constraints from SN light curves using Approximate Bayesian Computation.',
		license='MIT',

		classifiers=[
			'Development Status :: 5 - Production/Stable',
			'Environment :: Console',
			'Operating System :: OS Independent',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: MIT License',
			'Programming Language :: Python :: 2.7',
			'Topic :: Scientific/Engineering',
			],
		requires=['NumPy (>=2.7)',],
		long_description="""
		Approximate Bayesian computation (ABC) and so 
		called "likelihood free" Markov chain Monte Carlo 
		techniques are popular methods for tackling parameter 
		inference in scenarios where the likelihood is intractable or unknown. 
		These methods are called likelihood free as they are free from 
		the usual assumptions about the form of the likelihood e.g. Gaussian, 
		as ABC aims to simulate samples from the parameter posterior distribution directly.
		``superABC`` is a python package that implements  the astroABC sampler using either the SNANA or the sncosmo Supernovae simulation
		packages and simulates samples from the posterior distribution using two distinct distance metrics.
		``superABC`` requires ``NumPy``,``SciPy`` and ``astroabc``. The SNANA and rootpy packages are required if using this as a forward model simulation.
		sncosmo, astropy and pandas are required if using sncosmo as the forward model simulation. 
		""",

		#packages=find_packages(exclude=['contrib', 'docs', 'tests']),
		#packages=["astroab],
		)
		            
if __name__ == '__main__':
    setup_superabc()

