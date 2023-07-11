# functions to compute bulk information

#' Create functions that interpolate a spike profile, spike*a*b^2 + B
#'
#' @param spike numeric vector of spike values
#' @param singular_value spike singular values
#' @param u_proj2 square of projection of true singular vector (a) against estimated
bulk_profile_interpolation.analysis_util = function(spike, singular_value, u_proj2)
{
  spike_max = max(spike)

  singular_value_f = splinefun(x = spike,
                               y = singular_value)
  c_proj2_f = splinefun(x = spike,
                            y = u_proj2)

  # just in case there's an overflow
  singular_value_f_extended = function(a) {
    ifelse(a < spike_max, singular_value_f(a),
           a + singular_value_f(spike_max) - spike_max)
  }
  c_proj2_f_extended = function(a) {
    ifelse(a < spike_max, c_proj2_f(a), 1)
  }

  return (list(singular_value=singular_value_f_extended,
               coherent_proj2=c_proj2_f_extended))

}


#' Creates a bulk object, used in other util functions below
#'
#' @param s Seurat object
#' @param X the presumed mean expression
bulk.analysis_util = function(s, X)
{
  Y = scaled_counts.dataset(s)

  B = Y - X
  e = eigen(t(B) %*% B)
  V = e$vectors
  U = B %*% V
  U = apply(U, 2, function(x) x/sqrt(sum(x^2)))
  vals = e$values

  # remove the smallest eigenvalue/vectors since those are 0
  # and their computation is unstable
  active_vals = which(vals > 1E-4)
  U = U[,active_vals]
  V = V[,active_vals]
  vals = vals[active_vals]

  return (list(B=B, U=U, V=V, vals=vals))
}

#' Uses spectrode to compute the spectrum approximation of the bulk
#'
#' @param b a bulk object created from bulk.analysis_util
#' @param vanilla_MP should the spectrode object be construted assuming
#' all columns have the same norm, giving
#'
#' @return a spectrode object providing density and support of the MP spectrum
spectrode_MP.analysis_util = function(b, vanilla_MP=F)
{
  B = b$B
  if (vanilla_MP)
    par = param(1, w=1, gamma=ncol(B)/nrow(B))
  else
    par = param(colSums(B^2), w=1, gamma=ncol(B)/nrow(B))
  sp = spectrode(par)

  return (sp)
}

#' A convenience function for me to check bulk against MP assumption
show_bulk.analysis_util = function(s)
{
  Y = scaled_counts.dataset(s)
  Xs = 0*Y
  b = bulk.analysis_util(s, Xs)
  sp = spectrode_MP.analysis_util(b)
  evals = eigen(t(Y) %*% Y)$values
  plot.spectrode(sp, evals)
}

#' Computes the filtered statistics of spikes added to bulk
#'
#' @param b a bulk object, created by bulk.analysis_util
#' @param max_spike maximum spike strength (lambda in + lambda*u*v^T)
#' @param d_spike grid size in spike values to consider
spike_profile.analysis_util = function(b, spikes)
{
  lambda_values = spikes
  sp = spectrode_MP.analysis_util(b)
  s_stats = spike_statistics.spectrode(sp, lambda_values=lambda_values)
  return (list(spike=lambda_values,
               singular_value=s_stats$singular_values,
               u_proj2=s_stats$proj_u2,
               v_proj2=s_stats$proj_v2))
}


###################################################################
# snr

#' Computes single to noise ratio across a collection of genes
#'
#' @param s A Seurat object
#' @param cell_states a vector giving each cell one of two values,
#' reflecting the two cell states of the assumed module
#' @param genes the genes (names) in the module
#'
#' @details The data matrix of s is used to compute snr
#'
#' @return a named (genes) numeric vector of snr's
snr.analysis_util = function(s, cell_states, genes)
{
  if (length(unique(cell_states)) != 2)
    stop("I assume 2 cell states")

  states = cell_states %>% unique
  states1 = cell_states == states[1]
  states2 = cell_states == states[2]
  p1 = mean(cell_states == states[1])
  p2 = 1 - p1

  D = s@assays$RNA@data[genes,,drop=F]

  snr = apply(D, 1, function(g) {

    g1 = g[states1]; g2 = g[states2]
    dmu2 = (mean(g1) - mean(g2))^2
    V1 = var(g1); V2 = var(g2)
    snr = p1*p2*dmu2/(p1*V1 + p2*V2)

    return (snr)
  }) %>% setNames(genes)

  return (snr)
}

########################################################################
# filtering functions

filtering_measure.analysis_util = function(s, u_bar, spike_k, bulk_k)
{
  captured_A = captured_energy.analysis_util(s, u_bar, spike_k, bulk_k)
  A_captured_spike = captured_A$spike
  A_captured_total = captured_A$total
  A_captured_bulk = captured_A$bulk


  df = data.frame(A_spike=A_captured_spike,
                  A_bulk=A_captured_bulk,
                  A_total=A_captured_total,
                  stringsAsFactors = F)
  return (df)
}

captured_energy.analysis_util = function(s, u_bar, spike_k, bulk_k)
{
  uv = s@misc$uv
  u = uv$u; dv = uv$d

  u_proj = t(u) %*% u_bar
  A = (u_proj * dv)

  A_spike = (sum(A[1:spike_k]^2))
  A_total = (sum(A[1:bulk_k]^2))
  A_bulk = (A_total - A_spike)

  return (list(spike=A_spike, bulk=A_bulk, total=A_total))
}

#########################################
# functions to compute the svd of the signal added to the baseline modules

# I need to compute eigenvalues and eigenvectors for spike^2*u*u^T + \sum_i d_i^2 x_ix_i^T
# for fixed u, X, but varying spike values.

# Let W = spike^2*u*u^T + \sum_k d_k^2 x_k x_k^T
# Let M be the matrix for basis u,x_1,x_2,...,x_k (not orthonormal)
# In this basis it's easy to compute WW^T, call this matrix Z
# Let M=QR.  Then Q^TWW^TQ = RZR^-1.
# Compute Za = l*a.  Then b = Ra.  RZR^-1(b) = l*b.   Q^TWW^TQ(b) = l*b
# I can then compute WW^TQb = Qb, but I want u^TQb, which is b[1]

#' Creates an svd perturbation object
#'
#' @param u perturbation vector
#' @param d baseline singular values
#' @param x baseline singular vectors
#'
#' @details stores the Z and R matrices which are needed to compute
#' the projection of u onto the perturbed eigenvectors and the perturbed
#' eigenvalues
svd_perturbation.analysis_util = function(u, d, x)
{
  M = cbind(u, x)
  qr_out = qr(M)
  Q = qr.Q(qr_out);
  R = qr.R(qr_out)
  tMu = as.numeric(t(M) %*% u)

  # hardcode Z, only first row will vary depending on A
  rank = ncol(M)
  Z = matrix(0, nrow=rank, ncol=rank)
  diag(Z) = c(1,d^2)
  Z[1,] = tMu
  Z[,1] = tMu * c(1,d^2)

  return (list(R=R, Z=Z, Q=Q))
}

#' Computes the eigenvalues of spike^2uu^T + \sum_i d_i^2 x_ix_i^T and
#' the projection of u onto the eigenvectors
#'
#' @param spike spike strength of perturbation
#' @param sp a spike perturbation object
#'
#' @return a list containin eigenvalues and squared projection
compute_svd_perturbation.analysis_util = function(spike, sp)
{
  Q = sp$Q
  R = sp$R; Z = sp$Z
  Z[1,] = spike^2*Z[1,]
  e = eigen(Z)
  values = e$values
  vectors = R %*% e$vectors

  vectors = apply(vectors, 2, function(v) v/sqrt(sum(v^2)))
  u_dot = vectors[1,]

  full_eigenvectors = Q %*% vectors
  # debug, compare above against true values
  #W = A[i]*u%*% t(u) + Xs %*% t(Xs)
  #uv = irlba::svdr(W, k=rank+2)
  #W_values = uv$d
  #W_u_dot = t(uv$u) %*% u

  return (list(eigenvalues=values, u_proj2=u_dot^2,
               eigenvectors=full_eigenvectors))
}


#######################################################################
# clustering functions
# functions to computing theoretical clustering based on bulk and spike

#' Computes nn distribution for a spike that is not necessarily orthogonal to
#' the baseline
#' 
#' We now have a spike matrix where u is the signal and A is the signal strength.
#' Below lambda_i is the rho_i in the manuscript, the strength of a baseline module.
#'  X^TX = \sum lamdba_i^2 a_i a_i^T + A u u^T
#'  Applying a svd to the spike matrix
#'  X^TX = \sum rho_i^2 x_i x_i^T
#'  Passing these spikes through the bulk,
#'  \sum tilde(rho)_i^2 (alpha_i x_i + sqrt(1-alpha_i^2) noise) (...)^T
#'  Then we want to express in terms of u
#'  \sum tilde(rho)_i^2 (alpha_i beta_i u_i + alpha_i*sqrt(1-beta_i^2)tilde(x)_i
#'                      + sqrt(1-alpha_i^2) noise) (...)^T
#'  where tilde(x)_i is orthogonal to u.
#'  
#'  @param u The eigenvectors of the baseline PCA (not the signal)
#'  @param d The eigenvalues of the PCA
#'  @param celltype character vector giving perturbation state (one of "stim" and "ctrl")
#'  @param spike_sv_values the filtered signal strength, can be a vector for fast processing
#'  of many signal strenght values
#'  @param spike_k PCA dimension
#'  @param bulk_sv_f function computing singular value after passing through the bulk
#'  @param bulk_alpha2_f function computing alpha2 after passing through bulk.
nn_theory_none.analysis_util = function(u, d, celltype, spike_sv_values,
                                        spike_k,
                                        bulk_sv_f, bulk_alpha2_f,
                                        unlocalized,
                                        N=1,
                                        debug=0)
{
  if (spike_k==1) {
    message("code is written assuming spike_k > 1")
    stop(1)
  }

  ncell = length(celltype)
  stim_ind = which(celltype == "stim")
  ctrl_ind = which(celltype == "ctrl")
  p = sum(celltype =="stim")/ncell

  # form the module signal
  u_spike = make_u_spike.analysis.util(celltype)
  # apply svd to the module signal and the baseline spikes, sp
  # is an object that allows me to quickly compute the svd for 
  # different signal strengths
  sp = svd_perturbation.analysis_util(u_spike, d[1:spike_k,drop=F],
                                      u[,1:spike_k,drop=F])
  nn_vals = sapply(1:N, function(jj) {
   cat("simulating nn", jj, "of", N, "\n")
   sapply(spike_sv_values, function(c_sv) {

    if (debug == 1)
      cat("processing spike value", c_sv, "\n")
    if (debug == 2)
      browser()
     
    # compute the svd given a particular signal strength
    svd_info = compute_svd_perturbation.analysis_util(c_sv, sp)

    # new spike matrix,  \sum rho_i x_i y_i^T
    beta2_full = svd_info$u_proj2
    rho_full = sqrt(svd_info$eigenvalues)
    X_full = svd_info$eigenvectors

    # put the spikes in decreasing order and keep only spike_k
    end_ind = spike_k+1
    rho_ind = order(rho_full, decreasing=T)
    rho_full = rho_full[rho_ind]; beta2_full = beta2_full[rho_ind]
    rho = rho_full[1:end_ind]; beta2 = beta2_full[1:end_ind]
    X = X_full[,rho_ind]; X = X[,1:end_ind]

    # pass through the bulk
    tilde_rho = pmax(bulk_sv_f(rho),0)
    alpha2 = pmax(bulk_alpha2_f(rho),0)

    X_loc = sapply(1:end_ind, function(i) {
      a2 = alpha2[i]; b2 = beta2[i]
      # I need to recompute b2 to get the sign
      
      b2_s = sum(X[,i]*u_spike)
      #if (abs(abs(b2_s) - sqrt(b2)) > 1E-5)
      #  browser()

      v = sqrt(a2)*X[,i] + sqrt(1-a2)*rnorm(ncell)/sqrt(ncell)
      return (v)
    })

    inds = 1:end_ind
    nn_loc = nn_calculation.analysis_util(X_loc, tilde_rho, 
                                          inds, celltype, scale=ncell)

    if (debug==2)
      browser()

    return (nn_loc$nn)
   })
  })
  if (N > 1) {
    s2 = apply(nn_vals, 1, var)
    nn_vals1 = rowMeans(nn_vals)
  }
  else {
    nn_vals1 = as.numeric(nn_vals)
    s2 = 0
  }

  return (list(nn=nn_vals1, s2=s2))
}

#' Computes nn distributions for a spike assumed to be orthogonal
#' to a homogenous baseline
#' 
#' The function considers multiple perturbation spikes, but the
#' nn metric is computed for each using the homogenous baseline
#' assumption
#'
#' @param s Seurat object for the baseline dataset (assumed to svd in s@misc$uv)
#' @param spike_k PCA dimension
#' @param sv a vector of singular values (tilde(lambda) in the manuscript)
#' @param alpha2 vector of the squared projection of the spike eigenvector
#' onto the sampled eigenvector (alpha squared in the manuscript)
#' @param N number of times to sample from null for Monte Carlo estimate
#' 
#' @return a list containing the nearest neighbor metric estimate (nn)
#' and the variance of the estimator (s2)
nn_theory_homogeneous.analysis_util = function(s, spike_k, sv, alpha2,
                                          N=1,
                                          debug=F)
{
  if (spike_k==1) {
    message("code is written assuming spike_k > 1")
    stop(1)
  }

  d = s@misc$uv$d
  min_lambda = min(d[1:spike_k])

  ncell = ncol(s)
  n1 = sum(s$stim=="stim")
  n2 = sum(s$stim == "ctrl")
  stim_ind = which(s$stim=="stim")
  ctrl_ind = which(s$stim == "ctrl")

  p = n1/(n1+n2)

  nn = sapply(1:N, function(i) {
    mapply(function(c_sv, c_alpha2, i) {
      if (debug)
       cat("computing", c_sv, ": ", i, "of", length(sv), "\n")

      # although not needed, I split into whether alpha2 > 0
      
      # Case 1 : alpha2 > 0
      u0_theory = rep(0, ncell)
      u0_theory[stim_ind] = sqrt(c_alpha2)/sqrt(p*(1-p)) + rnorm(n1, sd=sqrt(1-c_alpha2))
      u0_theory[ctrl_ind] =  rnorm(n2, sd=sqrt(1-c_alpha2))
      u0_theory = u0_theory/sqrt(ncell)

      # Case 2 : alpha2 = 0
      u = matrix(rnorm(ncell*spike_k), nrow=ncell, ncol=spike_k)/sqrt(ncell)

      # Case 2, the spike is lost 
      if (c_sv < min_lambda)
          info = nn_calculation.analysis_util(u, d, 1:spike_k, s$stim, scale=ncell)
      # Case 1: the spike is in the PCA
      else
          info = nn_calculation.analysis_util(cbind(u0_theory, u),
                                          c(c_sv, d),
                                          1:spike_k,
                                          s$stim,
                                          scale=ncell)

      return (info$nn)

    }, sv, alpha2, 1:length(sv))
  })
  if (N > 1) {
    s2 = apply(nn, 1, var)
    nn = rowMeans(nn)
  } else {
    s2 = 0
  }

  return (list(nn=nn, s2=s2))
}



#' Calculate the nearest neigbhor metric from a Seurat object
#'
#' @param s Seurat object
#' @param spike_k number of pcs to use
#' @param N subsample size since it's too expensive to compute full distance matrix in general
nn_sampled.analysis_util = function(s, spike_k, N=3000)
{
  uv = s@misc$uv
  celltypes = s$stim
  ncell = length(celltypes)

  if (spike_k==1)
    D = uv$d[1:spike_k]*diag(1)
  else
    D = diag(uv$d[1:spike_k])
  W_full = uv$u[,1:spike_k,drop=F] %*% D

  N = min(N, nrow(W_full))
  active_ind = sample.int(nrow(W_full), N)

  W_full = W_full[active_ind,]
  stim = s$stim[active_ind]

  d_full = dist(W_full) %>% as.matrix
  diag(d_full) = Inf
  d_full = d_full^2

  w_nn = apply(d_full[stim=="stim", stim=="stim"], 1, min)
  a_nn = apply(d_full[stim=="stim", stim=="ctrl"], 1, min)

  return (mean(w_nn < a_nn))
}





#' Calculate the nearest neighbor metric from PCA information
#' 
#' @param u the PCA eigenvector matrix
#' @param lambda the PCA eigenvalues
#' @param celltype a vector of celltypes, containing "stim" and "ctrl"
#' @param N if the number of cells is very large, computing the exact
#' nearest neighbor metric is not feasible.  Instead, subsample N 
#' cells.
nn_calculation.analysis_util = function(u, lambda,
                                        active_ind,
                                        celltype,
                                        scale=NULL,
                                        N=3000)
{
  if (is.vector(u))
    u = matrix(u, ncol=1)

  u = u[,active_ind,drop=F]
  lambda = lambda[active_ind]
  if (length(lambda)==1)
    D = lambda*diag(1)
  else
    D = diag(lambda)
  W = u %*% D

  N = min(N, nrow(W))
  active_rows = sample.int(nrow(W), N)
  W = W[active_rows,]
  celltype = celltype[active_rows]

  d = dist(W) %>% as.matrix()
  diag(d) = Inf
  d = d^2
  if (!is.null(scale))
    d = scale*d

  stim_ind = which(celltype == "stim")
  ctrl_ind = which(celltype == "ctrl")

  w_nn = apply(d[stim_ind, stim_ind], 1, min)
  a_nn = apply(d[stim_ind, ctrl_ind], 1, min)

  nn = mean(w_nn < a_nn)

  return (list(w_nn=w_nn, a_nn=a_nn, d=d,
               nn=nn))
}

make_u_spike.analysis.util = function(celltype)
{
  ncell = length(celltype)
  v = rep(0, ncell)
  v[celltype=="stim"] = 1
  v = v - mean(v)
  v = v/sqrt(sum(v^2))
  return (v)
}



