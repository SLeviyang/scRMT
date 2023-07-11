#' Compute the s(z) the Stieljes transform, over a range of real z that are
#' outside the support intervals of the eigenvalue density.  Useful
#' in computation of spike critical points and other statistics.
#'
#' @param z0 a scalar, start of interval.  Must be outside the support
#' @param z1 a scalar, end of interval.  Must be outside the support
#' @param support.  A k x 2 matrix giving intervals of support of ESD.  This
#' can be the support field of a spectrode object, but the density computed
#' in the construction of a spectrode object is not needed.
#' @param param a spectrode paramater object
#' @param dz grid width
#'
#' @return a vector of imaginary values giving the Stieljes value on [z0,z1]
Stieljes = function(z0, z1, support, param, dz=1E-5)
{
  if (z1 < z0) {
    message("z1 must be greater than z0")
    return (NULL)
  }
  # derivative of s(z)
  mp = function(z,m,parms) {
    t = parms$t; w = parms$w; g = parms$gamma
    denom = 1/m^2 - g*sum(w*t^2/(1+m*t)^2)
    return (list((z0-z1)/denom))
  }

  # make sure the interval is outside the support
  outside = apply(support, 1, function(interval) {
    a = interval[1]; b = interval[2]
    if (z1 < a | z0 > b)
      return (T)
    return (F)
  })

  if (!all(outside)) {
    message("[z0,z1] must be outside the support of the bulk density")
    return (NULL)
  }

  #find starting point for ODE using hybrid fixed point method:
  zh = z1 + 10*1i

  m1 = compute_Stieltjes_hybrid(zh, z1, param, maxIter = 1E6, dz=1E-5)
  grid = seq(0, 1, dz)
  ode_out = deSolve::zvode(m1, times=grid, func=mp, parms=param,
                             rtol=param$tol, atol=param$tol)

  return (list(grid=grid*z0 + (1-grid)*z1, m=Re(ode_out[,2])))
}

