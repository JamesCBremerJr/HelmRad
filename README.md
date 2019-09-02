A fast solver for the variable coefficient Helmholtz equation in the radially symmetric case
=============================================================================================

Standard solvers for the variable coefficient Helmholtz equation in two spatial
dimensions have running times which grow quadratically with the wavenumber k
This package contains a solver whose complexity is quasilinear in the wavenumber k,
but which only applies the case of a radially symmetric scattering potential.
