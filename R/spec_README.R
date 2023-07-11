# files

# spec_spectrode : constructor functions for param and spectrode objects.
#. Construction of spectode depends on support.R and density.R

#' spec_support : functions to compute the support
#' spec_density : once support is computed, computes the density using ODE
#' spec_fixed_point : fixed point computations used as starting point of ODE.
#' Include a hybrid function in which a fixed point is computed and then ODE
#' is integrated to get to another point, without returning grid.
