#' Construct a spectrodeR parameter object
#'
#' @param t a numeric vector giving the population covariance
#' eigenvalues
#' @param w a numeric vector giving weights for the t variable
#' entries.  Automatically normalized to 1.  w can be a single value,
#' in which case it is recycled, otherwise it must have the same length as t.
#' @param gamma the matrix dimension ratio, i.e. p/n, ncols/nrows
#' @param tol tolerance at which computations, mostly root finding,
#' are executed
#' @param merge_tol tolerance for differences in values of t.
#'
#' @return a spectrodeR object.
#' @export
param = function(t, w, gamma, tol=1E-6,
                 merge_tol=1E-3)
{
  if (length(w)==1)
    w = rep(1, length(t))
  if (length(w) != length(t)) {
    message("w must be length 1 or the length of t")
    return (NULL)
  }
  if (any(w <= 0) | any(t < 0)) {
    message("all t must be nonnegative and w must be positive")
    return (NULL)
  }
  w = w/sum(w)

  ord_ind = order(t)
  t = t[ord_ind]
  w = w[ord_ind]

  # iteratively merge any value of t within merge_tol
  nt = length(t)
  new_t = NULL; new_w = NULL
  base_ind = 1; cur_ind = 2
  while(base_ind <= nt) {
    while(t[base_ind] > (t[cur_ind] - merge_tol) & cur_ind <= nt)
      cur_ind = cur_ind + 1

    new_t = c(new_t, mean(t[base_ind:(cur_ind-1)]))
    new_w = c(new_w, sum(w[base_ind:(cur_ind-1)]))
    base_ind = cur_ind
  }

  p = list(t=new_t, w=new_w, gamma=gamma, tol=tol)
  class(p) = c("sprectrodeR parameter", "list")

  return (p)
}


#' Compute spectral density of the sample covariance matrix in the
#' large n,p limit.
#'
#' @param param a spectrodeR parameter object
#' @param dx grid width for the density
#' @param epsilon steiltjes transform are computed at (x + i*epsilon)
#'  where x is in the support.  Smaller episilon are more accurate
#'  for computing ESD but require more iterations.
#'  @param plot should the density be plotted
#'
#' @return a list containing fields
#'  @field grid grid of x values
#'  @field density density (ESD) at grid points
#'  @field support a K by 2 matrix containing the K support intervals
#'  @filed dx the grid width (convenient for integrating against the density)
#'
#' @export
spectrode = function(param, dx=1E-4, epsilon=1E-6)
{
  I = support(param)
  d = density(I, dx, param, epsilon=epsilon)

  return (list(grid=d$grid,
               density=d$density,
               support=I,
               dx=dx,
               param=param))
}
