#' These functions consider a the spike model of
#'  Benaych-Georges and Nadakuditi 2012, JMV
#'
#' Let X be an n by N random matrix with the N by N sample covariance matrix
#' spectrum (ESD) given by the spectrodeR param density.
#' We consider the spiked model,
#'  A =  X + lambda*u*v^T
#' where u and v are vectors of dimension n and N, respectively.

#' Computes the critical point at which a spike will exit the bulk
#'
#'  Results by Benaych-Georges and Nadakuditi, and many other authors,
#'  show that the dominant eigenvalue of A^TA will be the same as X^TX for
#'   lambda^2 less than a detection threshold, computed by this function.
#'
#' @param s a spectrode object
#' @param epsilon a regularization constan
#'
#' @return critical value for lambda
#'
#' @details This function computes the critical points using the formula,
#' in Benaych-Georges and Nadakuditi 2012, JMV
#'

#' detection threshold = 1/sqrt(D(b^2)), where
#'
#' D(z) = \int_S z/(z - l) dF_X^TX(l) * \int_S z/(z - l) dF_X^TX(l)
#' where dF gives ESD, as provided by the spectrode object and
#' where b^2 is the outer edge of support.
#'
#' See Leeb 2021 Adv Comput Math for an
#' alternative method of computing the detection threshold
#' that doesn't require the ESD.
spike_detection_threshold.spectrode = function(s, epsilon=1E-5)
{

  b2 = max(s$support) + 1E-5
  z_vals = sqrt(b2)
  f = compute_spike_functions.spectrode(s, z_vals, return_d = F)
  D = f$D

  return (c(threshold=1/sqrt(D)))
}

#' Computes projection of dominant eigevalue onto the spike vectors u,v over
#' a grid of lambda values
#'
#' @param spectrode object
#' @param lambda_values values of lambda for which statistics are to be computed
spike_statistics.spectrode = function(s, lambda_values)
{
  lambda_values = sort(lambda_values)
  threshold = spike_detection_threshold.spectrode(s)

  # remove anything less than threshold, and add trivial stats at the end
  pre_threshold = lambda_values <= threshold
  lambdas = lambda_values[!pre_threshold]
  if (length(lambdas > 0)) {
    # find the singular value for each lambda.
    # I need to solve D(z) = 1/lambda^2
    # create a grid
    b2 = max(s$support) + 1E-10
    z_start = sqrt(b2)
    z_end = z_start+1
    while (compute_spike_functions.spectrode(s, z_end)$D > min(1/lambdas^2))
      z_end = z_end + 2

    # I need a lot of resolution near the critical point
    z_vals = c(seq(z_start, z_start+.2, 0.0001), seq(z_start+.1+0.01, z_end, 0.001))
    print(length(z_vals))
    D_vals = compute_spike_functions.spectrode(s, z_vals, return_d=F)$D

    # solve D(z) = 1/lambda^2
    sv = sapply(lambdas, function(l) {
      ind = which.min(abs(D_vals - 1/l^2))
      z_vals[ind]
    })

    f_vals = compute_spike_functions.spectrode(s, sv, return_d=T)
    phi_XTX = f_vals$phi_XTX
    phi_XXT = f_vals$phi_XXT
    d_D = f_vals$d_D

    proj_u2 = -2*phi_XXT/lambdas^2/d_D
    proj_v2 = -2*phi_XTX/lambdas^2/d_D
  }

  sv_all = proj_u2_all = proj_v2_all = rep(NA, length(lambda_values))
  sv_all[pre_threshold] = 0
  proj_v2_all[pre_threshold] = 0
  proj_u2_all[pre_threshold] = 0

  sv_all[!pre_threshold] = sv
  proj_v2_all[!pre_threshold] = proj_v2
  proj_u2_all[!pre_threshold] = proj_u2


  return (list(lambda=lambda_values,
               singular_values=sv_all,
               proj_u2=proj_u2_all, proj_v2=proj_v2_all))
}


#' Compute the D-transform and related functions
#'
#' @param z_vals point at which to evaluate functions
#' @param return_d should the derivative functions also be computed
compute_spike_functions.spectrode = function(s, z_vals,
                                             return_d=F)
{
  param = s$param
  gamma = param$gamma
  density = s$density
  x = s$grid
  dx = s$dx

  total = sum(density*dx)
  density = density/total
  if (gamma > 1)
    g = 1/gamma
  else
    g = gamma

  phi_XTX = sapply(z_vals, function(z) sum(z/(z^2 - x)*density*dx))
  phi_XXT = sapply(z_vals, function(z)
              sum(g*z/(z^2 - x)*density*dx) + (1-g)/z)
  D = phi_XTX*phi_XXT


  if (return_d) {
    d_phi_XTX = sapply(z_vals, function(z)
                sum((1/(z^2-x) - 2*z^2/(z^2 - x)^2)*density*dx))
   d_phi_XXT = sapply(z_vals, function(z) {
    sum(g*(1/(z^2-x) - 2*z^2/(z^2 - x)^2)*density*dx) - (1-g)/z^2
    })
  }

  if (return_d)
    d_D = d_phi_XTX*phi_XXT + phi_XTX*d_phi_XXT
  else
    d_D = NULL

  return (list(phi_XTX=phi_XTX, phi_XXT=phi_XXT, D=D, d_D=d_D))
}
