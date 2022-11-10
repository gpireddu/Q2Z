Q2Z
==========

Q2Z.py is a Python code for computing the frequency-dependent admittance / impedance spectra from the time series of the total electrode charges.

This method is based on linear response theory and relies on the Fourier-Laplace transform of the time autocorrelation function of the total charges fluctuations.

---
# Reference

[Giovanni Pireddu and Benjamin Rotenberg,arXiv:2206.13322](https://doi.org/10.48550/arXiv.2206.13322)

---
# Usage
The electrode charge time series should be in the same folder and named ```total_charges.out```. An example is provided in ```Example/total_charges.out``` for testing purposes only.

Simply run:

```python3 Q2Z.py```

---
# Notes
* The code assumes a MetalWalls output format [https://gitlab.com/ampere2/metalwalls/-/wikis/output-files#total_charges.out](https://gitlab.com/ampere2/metalwalls/-/wikis/output-files#total_charges.out)
* Parameters such as the temperature, duration of time steps and the list of desired frequencies are hard-coded, for the moment.
* The code is not tested for Python versions other than Python3

---
# Dependencies
Q2Z requires the following packages to run:
* numpy
* scipy
* numba
