# TIPTOP

In order to be able to easily predict the AO performance, we have developed this fast algorithm producing the expected Adaptive Optics (AO see https://en.wikipedia.org/wiki/Adaptive_optics) Point Spread Function (PSF) for any of the existing AO observing modes (Single-Conjugate-AO, Laser-Tomographic-AO, Multi-Conjugate-AO, Ground-Layer-AO), and any atmospheric conditions. This TIPTOP tool takes its roots in an analytical approach, where the simulations are done in the Fourier domain. This allows to reach a very fast computation time (few seconds per PSF with GPU acceleration), and efficiently explore the wide parameter space.

A documentation under development can be found here https://tiptopdoc.readthedocs.io/en/latest/

## References

Reference: "TIPTOP: a new tool to efficiently predict your favorite AO PSF" SPIE 2020 (ARXIV: https://doi.org/10.48550/arXiv.2101.06486).

## Installation

### Anaconda setup

We advise to use conda.

Create a Anaconda environment for TIPTOP with python 3.X (below for example 3.8):
```
conda create -–name tiptop python=3.8

conda activate tiptop
```
Install cupy (it requires a CUDA-Enabled NVIDIA GPU):
```
conda install -c conda-forge cupy
```
Install other libraries:
```
conda install ipython matplotlib scipy astropy sympy
```
and to get jupiter and jupiter-lab (many files are Jupyter Notebooks):
```
conda install jupyter
conda install -c conda-forge jupyterlab
```

### Get the library

There are two ways clone TIPTOP:

1. with sub-modules:
```
git clone --recurse-submodules https://github.com/astro-tiptop/TIPTOP.git
```
2. without sub-modules (cloned separately):
```
git clone https://github.com/astro-tiptop/TIPTOP.git
```

In the second case please clone also the repos:

- https://github.com/astro-tiptop/MASTSEL
- https://github.com/astro-tiptop/SYMAO
- https://github.com/astro-tiptop/SEEING
- https://github.com/astro-tiptop/P3/

as reported in the `setup.py` file.

### Install with pip

Navigate to the folder in which you have cloned TIPTOP, then navigate to the
folder P3 and install P3 (remove `--user` to install to all users):
```
pip install --user -e .
```
Navigate out of P3, navigate inside MASTSEL and install MASTSEL:
```
pip install --user -e .
```
Navigate out of MASTSEL, navigate inside SEEING and install SEEING:
```
pip install --user -e .
```
Navigate out of SEEING, navigate inside of SYMAO and install SYMAO:
```
pip install --user -e .
```
Navigate out of SEEING and install TIPTOP:
```
pip install --user -e .
```
