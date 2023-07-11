#' Compute support of spectrum
#'
#' @param param a spectrodeR parameter object
#'
#' @returns a K by 2 matrix containing the support intervals
#' @export
support = function(param)
{
  c_intervals = support_C(param) %>%
                dplyr::select(z1, z2) %>% as.matrix
  K = nrow(c_intervals)
  intervals = matrix(0, nrow=K-1, ncol=2)
  for (i in 1:(K-1))
    intervals[i,] = c(c_intervals[i,2], c_intervals[i+1,1])

  return (intervals)
}

#' Calculate z(s), z'(s), z''(s), ....
#'
#' @param s a vector of s values
#' @param k index of the derivative, i.e. z(s) is k=0, z'(s) is k=1
#' @param param t,w,gamma, see help for compute_spectrum_support
#'
#' @returns a numeric vector
z_func = function(s, k, param) {
  t = param$t; w = param$w; g = param$gamma

  t_s = t %*% t(s)
  mat = 1/(1 + t_s) * t

  kp1 = k + 1
  term = -1/s^kp1 + g*colSums(mat^kp1*w)
  return (factorial(k)*(-1)^k*term)
}

#' Determine the complement of the support of the spectrum
#'
#' @param param a sprectrodeR parameter object
#'
#' @returns a data.frame containing the s (stieljes transform) and
#' z values that specify the intervals for the complement of the
#' support
#'
#' @details Uses the ideas mentioned in the original Marchenko-Pastur paper
#' along with the further details in Silverstein, Choi 95.  Essentially,
#' the complement is contained in z'(s) > 0.  The function computes the critical
#' points z'(s) = 0 which allows for the construction of the support.
#' @export
support_C = function(param)
{
  t = param$t; gamma = param$gamma
  # First, collect the edge cases, to the left and right of -1/tmin, -1/tmax.
  # These cases differ according to gamma <= 1.
  if (gamma <= 1)
     intervals = support_C_edge_gamma_small(param)
  else
     intervals = support_C_edge_gamma_big(param)

  if (length(t)==1) {
    intervals = dplyr::arrange(intervals, z1)
    return (intervals)
  }

  # Now collect all supports between the -1/t_i
  t_sort = sort(t); n_intervals = length(t) - 1
  z_intervals = matrix(0, nrow=0, ncol=2)
  for (i in 1:n_intervals) {
    new_interval = support_C_between_ti(t_sort[i], t_sort[i+1], param)
    if (!is.null(new_interval))
      intervals = rbind(intervals, new_interval)
  }

  # put intervals in order
  intervals = dplyr::arrange(intervals, z1)

  return (intervals)
}

#' Compute the complement support within [-1/t_i, -1/t_{i+1}]
#'
#' @param t1,t2 points in the spectrum of H.  Must have t1 < t2.
#' @param p SpectrodeR parameter object
#'
#' @details  On [-1/t_i,-1/t_{i+1}], z''(s) = 0 has a unique solution,
#   implying that z'(s) = 0 has two solutions in the interval or none,
#   see Silverstien, Choi's Theorem 4.3.
support_C_between_ti = function(t1, t2, p)
{
  tol = p$tol

  s1 = -1/t1; s2 = -1/t2
  # if t1 and t2 are very close, we won't have the machine precision
  if(abs(s2-s1) < 1E-10)
    stop("the population eigenvalues are too close!")
  ds = abs(s2-s1)/100
  s1 = s1 + ds
  s2 = s2 - ds

  # find z''(s_c) = 0
  sc = uniroot(z_func, interval=c(s1,s2), k=2,
               p, tol=tol)$root
  # two z'(s) = 0 roots if z'(sc) > 0, otherwise no roots
  if (z_func(sc, k=1, p) > 0) {
    left_root = uniroot(z_func, interval=c(s1,sc), k=1,p,tol=tol)$root
    right_root = uniroot(z_func, interval=c(sc,s2), k=1,p,tol=tol)$root
    return (z_s_interval(left_root, right_root, p, "ti"))
  } else
    return (NULL)
}

#' Find the edge cases of for the support of the complement spectrum when gamma <= 1
#'
#' @param t,w,gamma,tol see help for compute_spectrum_support
#'
#' @returns a 2 by 2 vector containing the left and right most portions
#' of the support complement as a function of s.
#'
#' @details The two intervals [-1/t_max, 0] and [-Inf, -1/t_min] need
#' to be considered.
support_C_edge_gamma_small = function(param)
{
  tol = param$tol
  t = param$t

  # First the [-Inf, -1/t_min] support
  tmin = min(t)
  sR = -1/tmin - tol
  sL = -1/tmin - 1
  while (z_func(sL, k=1, param) < 0)
    sL = 1.2 * sL
  sc = uniroot(z_func, interval=c(sL, sR), k=1, param, tol=tol)$root
  intervals = z_s_interval(-Inf, sc, param, "edge1")

  # interval within [-1/tmax,0]
  sL = -1/max(t) + tol
  sR = -tol
  sc = uniroot(z_func, interval=c(sL, sR), k=1, param, tol=tol)$root
  intervals = rbind(intervals,
                    z_s_interval(sc, 0, param, "edge2"))

  return (intervals)
}




#' Find the edge cases of for the support of the complement spectrum when gamma > 1
#'
#' @param t,w,gamma,tol see help for compute_spectrum_support
#'
#' @returns a 2 by 2 vector containing the left and right most portions
#' of the support complement as a function of s.
#'
#' @details The two intervals [-1/t_max, 0] and [0,Inf] need
#' to be considered.
support_C_edge_gamma_big = function(param)
{
  t = param$t; tol = param$tol

  # interval within [-1/tmax,0]
  sL = -1/min(t) + tol
  sR = -tol
  sc = uniroot(z_func, interval=c(sL, sR), k=1, param, tol=tol)$root
  intervals = z_s_interval(sc, 0, param, "edge1")

  # Need to bound the [0,Inf] interval before root finding
  sL = tol
  sR = 1
  while (z_func(sR, k=1, param) > 0)
    sR = 1.2 * sR
  sc = uniroot(z_func, interval=c(sL, sR), k=1, param, tol=tol)$root
  intervals = rbind(intervals,
                    z_s_interval(0, sc, param, "edge2"))

  return (intervals)
}

#' Converts an s interval to a z interval
#'
#' @param s1,s2 end points of s intervals
#' @param p a param object
#' @param label a string
#'
#' @return a data.frame with one row containing fields s1, s2, z1, z2, label
#'
#' @details if s==0, then z is either Inf or -Inf depending on whether
#' the interval is negative/positive
z_s_interval = function(s1, s2, p, label=NA)
{
  if (s1==0)
    z1 = -Inf
  else
    z1 = z_func(s1, k=0, p)
  if (s2 == 0)
    z2 = Inf
  else
    z2 = z_func(s2, k=0, p)

  return (data.frame(s1=s1, s2=s2,
                     z1=z1, z2=z2,
                     type=label))
}

