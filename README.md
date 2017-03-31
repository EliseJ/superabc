<a href="https://github.com/EliseJ/superabc"><img src="https://github.com/EliseJ/superabc/blob/master/superabc_logo.001.jpeg"
align="left" hspace="5" vspace="3"></a>





<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br> <br> <br><br><br><br><br><br><br><br><br><br><br>

<br>

[![Latest Version](http://img.shields.io/pypi/v/superabc.svg?style=flat)](https://pypi.python.org/pypi/superabc/)
[![Open Source Love](https://badges.frapsoft.com/os/mit/mit.svg?v=102)](https://github.com/EliseJ/superabc/blob/master/LICENSE.txt)
 [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/EliseJ/superabc/issues)



Author: Elise Jennings

[arXiv:1611.03087](https://arxiv.org/abs/1611.03087)


superABC is a new sampling method for obtaining cosmological constraints from SN light curves using Approximate Bayesian Computation (ABC).

superABC comes with an interface to two forward model simulations for SN cosmology
* [sncosmo](http://sncosmo.readthedocs.io/en/v1.3.x/)
* SNANA

Please email elise@fnal.gov if you have any questions. 


### Installing ###

Install superABC using pip

```
$ pip install superabc
```

or git clone the repository using the url above. 
Check the dependencies listed in the next section are installed.

### Dependencies ###

* numpy
* scipy
* astroabc
* mpi4py
* multiprocessing
* sklearn

if using sncosmo:

* sncosmo
* pandas
* astropy

if using SNANA:

* SNANA
* rootpy


### License ###

Copyright 2016 Elise Jennings

superABC is free software made available under the MIT License. For details see the LICENSE.txt file.
