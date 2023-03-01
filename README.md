Q2Z
==========

Q2Z.py is a Python code for computing the frequency-dependent admittance / impedance spectra from the time series of the total electrode charges.

This method is based on linear response theory and relies on the Fourier-Laplace transform of the time autocorrelation function of the total charges fluctuations.

---
# Reference

[Giovanni Pireddu and Benjamin Rotenberg,Physical Review Letters, 130, 098001, 2023](https://doi.org/10.48550/arXiv.2206.13322)

Bibtex:
```
@article{PhysRevLett.130.098001,
  title = {Frequency-Dependent Impedance of Nanocapacitors from Electrode Charge Fluctuations as a Probe of Electrolyte Dynamics},
  author = {Pireddu, Giovanni and Rotenberg, Benjamin},
  journal = {Phys. Rev. Lett.},
  volume = {130},
  issue = {9},
  pages = {098001},
  numpages = {7},
  year = {2023},
  month = {Feb},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.130.098001},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.130.098001}
}
```

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
