Q2Z
==========

Q2Z.py is a Python code for computing the admittance / impedance spectra from the time series of the total electrode charges.

# Reference

[G. Pireddu and Benjamin Rotenberg,arXiv:2206.13322](https://doi.org/10.48550/arXiv.2206.13322)

# Notes
* The code assumes a MetalWalls output format [https://gitlab.com/ampere2/metalwalls/-/wikis/output-files#total_charges.out](https://gitlab.com/ampere2/metalwalls/-/wikis/output-files#total_charges.out)
* Parameters such as the temperature, duration of time steps and the list of desired frequencies are hard-coded, for the moment.
* The code is not tested for Python versions other than Python3

# Dependencies
Q2Z requires the following packages to run:
* numpy
* scipy
* numba
