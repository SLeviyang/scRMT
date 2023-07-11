#' Compute Stieljes transform at a point
#'
#' @param t population eigenvalues as non-negative numeric vector
#' @param w weight for the eigenvalues.  For example (t=1, w=1) gives MP distribution.
#' @param gamma aspect ratio, i.e. p/n.  Must be less than 1.
#' @param grid points at which spectrum should be evaluated
#' @param epsilon imaginary value to be added to the grid points
#'
#' @return a list with following fields
#'  @field density computed density at grid points
#'  @field m Stieljes transform on the grid
#'  @field numIter number of iterations performed
#'  @field stepSize
#'  @field grid the grid passed into the function, or the default created
#'
#'  @details
compute_Stieltjes_fixed_point <- function(z, param, maxIter=1E5)
{
  if (Mod(z)==0)
    z = 1E-8*1i
  t = param$t; w = param$w; gamma = param$gamma
  tol = param$tol

  Silverstein_func = function(m,z) -1/(z - gamma*sum(w*t/(1+t*m)))

  current_m <- Silverstein_func(-1/z,z)
  iter=0
  prev_m = Inf
  while(Mod(current_m - prev_m) > tol & iter < maxIter) {
    prev_m = current_m
    current_m = Silverstein_func(current_m, z)
    iter = iter + 1
  }

  return (list(m=current_m, iter=iter,
               error=Mod(current_m - prev_m)))
}


#' Compute m(z1) by using fixed point at m(z0) and then integrating to z1.
compute_Stieltjes_hybrid <- function(z0, z1, param, maxIter=1E5,
                                     dz=1E-5)
{
  # make sure z's are complex type
  z0 = z0 + 1i*0; z1 = z1 + 1i*0
  m0 = compute_Stieltjes_fixed_point(z0, param, maxIter=maxIter)$m

  mp = function(z,m,parms) {
    t = parms$t; w = parms$w; g = parms$gamma
    denom = 1/m^2 - g*sum(w*t^2/(1+m*t)^2)
    return (list((z1-z0)/denom))
  }

  grid = seq(0,1,dz)
  ode_out = deSolve::zvode(m0, times=grid, func=mp, parms=param,
                           rtol=param$tol, atol=param$tol)

  return (ode_out[nrow(ode_out),2])
}


