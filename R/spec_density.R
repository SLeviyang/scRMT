#' Compute the limit empirical spectral distribution given support
#'
#' @param support a K by 2 matrix containing the K intervals
#' @param dx grid width for density within each support interval
#' @param param a spectrodeR parameter object
#' @param epsilon imaginary portion of z
#'
#' @returns
#'  @field grid
#'  @field density
#' @export
density <- function(support, dx, param,
                      epsilon = 1E-6)
{
  # derivative of s(z)
  mp = function(z,m,parms) {
    t = parms$t; w = parms$w; g = parms$gamma
    denom = 1/m^2 - g*sum(w*t^2/(1+m*t)^2)
    return (list(1/denom))
  }

  # make sure support is ordered
  support = support[order(support[,1]),,drop=F]

  interval_info = apply(support, 1, function(cur_int) {
   #find starting point for ODE using hybrid fixed point method:
   z1_real = cur_int[1] + dx
   zh = z1_real + 10*1i
   z1 = z1_real + 1i*epsilon

   m1 = compute_Stieltjes_hybrid(zh, z1, param, maxIter = 1E7)
   grid = seq(z1_real, cur_int[2], by=dx)
   ode_out = deSolve::zvode(m1, times=grid, func=mp, parms=param,
                            rtol=param$tol, atol=param$tol)

   density = 1/pi*Im(ode_out[,2])
   return (list(grid=grid, density=density))
  })

  grid = sapply(interval_info, "[[", "grid") %>% unlist
  density = sapply(interval_info, "[[", "density") %>% unlist

  return (list(grid=grid, density=density))

}
