# Given full data and spike, compute information relating to the bulk
bulk.spike = function(s, X=NULL)
{
  if (is.null(X))
    X = s@misc$Xs %>% t
  Y = scaled_counts.dataset(s)

  G = Y - X
  p = param(colSums(G^2), w=1, gamma=ncol(X)/nrow(X))
  # compute the density for GG^T (n by n)
  spec = spectrode(p, dx=0.001)

  density = spec$density
  grid = spec$grid
  g = p$gamma
  dx = spec$dx

  if (g > 1) {
    message("I still need to code this case")
  }

  plot(grid, density)
  # I can compute the spike strength needed to exceed the bulk using
  D_info = D(max(spec$support) + 1E-6, max(spec$support)+4, spec$support, p)
  A_transition = 1/max(D_info$D)
  plot(D_info$grid, D_info$D)

  # solve for spike strength that leads to sv above the bulk
  return (A_transition)
}

# For each gene module, compute the cell singular vector using the spike matrix
module_eigenvectors.spike = function(s)
{
  Xs = s@misc$Xs %>% t
  gp = s@misc$gene_patterns %>% unique
  m = sapply(gp, function(pattern) {
    Z = Xs[,s@misc$gene_patterns == pattern]
    uv = irlba::svdr(Z, k=1)
    return (uv$d*uv$u)
  })
  colnames(m) = gp
  return (m)
}


#############################################

#' Construct the spike matrix after Seurat normalization from
#' the raw mean expression Matrix
Seurat_spike = function(X, dispersion, L=10000,
                        lib_sizes=NULL,
                        debug=F)
{
  # Seurat will normalize the genes (rows) of Y, so we need a
  # normalized version of X.  The normalization requires
  # us to know the variance of the noise, E[(Y_ij - X_ij)^2].

  # Seurat will consider log(1+max(Yij)/Lj*L)
  # We need 1) E[log(1 + NB(Xij,theta)/Lj*L)]
  #.        2) V[log(1 + NB(Xij,theta)/Lj*L)]
  pp = read_spike_processing(dispersion)

  logX_mean = matrix(0, nrow=nrow(X), ncol=ncol(X))
  logX_var = matrix(0, nrow=nrow(X), ncol=ncol(X))

  # I use the mean libsizes. Seurat will include noise in the libsizes
  # since it computes library size from Y
  if (!is.null(lib_sizes))
    lib_ratios = lib_sizes/L
  else
    lib_ratios = colSums(X)/L

  if (debug)
    browser()

  lib_ratio_bins = .bincode(lib_ratios, pp$lib_ratios,
                            right=F, include.lowest = T)
  X_bins = matrix(.bincode(as.numeric(X), pp$raw_expr_means,
                           right=F, include.lowest = T),
                  nrow=nrow(X), ncol=ncol(X))

  for (j in 1:ncol(X)) {
    cat("processing cell", j, "of", ncol(X), "\n")
    lr = lib_ratios[j]
    lr_index = lib_ratio_bins[j]
    if (is.na(lr_index)) {
      cat("library ratio outside of preprocessed range", lr, "\n")
    }

    for (i in 1:nrow(X)) {
      xij = X[i,j]
      if (xij == 0)
        next

      xij_index = X_bins[i,j]
      if (!is.na(xij_index) & !is.na(lr_index)) {
        logX_mean[i,j] = pp$logX_mean[xij_index, lr_index]
        logX_var[i,j] = pp$logX_var[xij_index, lr_index]
      } else {
        if (!is.na(lr_index)) {
          cat("expression outside of preprocessed range", xij, "\n")
        }
        Xstats = compute_logX_distribution(xij, dispersion, lr)
        logX_mean[i,j] = Xstats$logX_mean
        logX_var[i,j] = Xstats$logX_var
      }
    }
  }


  # Let Yi be expression of gene i.
  # Normalization :  (Yi - mean(Xi))/sqrt(\sum_j E[(Yij - Xi)^2])
  #        = (Yi - mean(Xi))/sqrt(sum_j E[(Y_ij - X_ij + Xij - Xi)])
  #.       = (Yi - mean(Xi))/sqrt(sum_j V[Y_ij] + (Xij - Xi)^2)
  # Form X normalized, keep gene_means and gene_normalization for analysis

  gene_means = rowMeans(logX_mean)
  gene_normalizations = sapply(1:nrow(X), function(i) {
    x = logX_mean[i,]
    mu = gene_means[i]
    sqrt(sum(logX_var[i,] + (x-mu)^2))
  })

  Xn = sapply(1:nrow(X), function(i) {
    if (gene_normalizations[i] > 0)
      (logX_mean[i,] - gene_means[i])/gene_normalizations[i]
    else
      rep(0,ncol(X))
  }) %>% t

  rownames(Xn) = rownames(X)
  colnames(Xn) = colnames(X)

  return (Xn)
  # return (list(Xs=Xn, gene_means = gene_means,
  #              logX_mean=logX_mean,
  #              gene_normalizations=gene_normalizations))
}


#' @param lib_ratio_parameters Gives a range of possibilities for Lj/L where Lj is
#' the library size for a cell and L is the Seurat normalization, typically
#' 10000.  Should be a numeric vector of length 3 giving min, max, and
#' number of grid points.
#' @param raw_expr_means_parameters A range of possible raw (i.e. not log)
#' expression values.  Should be a numeric vector of length 3 giving
#' min, max, and number of grid points.
#' @param dispersion the dispersion value for the NB model
#'
#' Seurat will normalize to log(1+max(Yij)/Lj*L) where Yij=NB(Xij,theta)
#' We need 1) E[log(1 + NB(Xij,theta)/Lj*L)]
#'         2) V[log(1 + NB(Xij,theta)/Lj*L)]
#' in order to compute the spike
create_spike_preprocessing = function(lib_ratios=c(seq(1/100,2,length.out=200),
                                                   seq(2.1,40,.1)),
                                      raw_expr_means=c(seq(0,1,.001),
                                                       seq(1.001,5,.001),
                                                       seq(5.01,300,.1),
                                                       seq(301,1000,2),
                                                       seq(1001, 10000, 20)),
                                      dispersion)
{
  sdir = file.path("spike_preprocessing")
  if (!dir.exists(sdir))
    dir.create(sdir)

  rds_file = file.path(sdir,
                       paste("spike_pre_dispersion_",  dispersion, ".rds",
                             sep=""))

  nlr = length(lib_ratios)
  nm = length(raw_expr_means)

  logX_mean = matrix(0, nrow=nm, ncol=nlr)
  logX_var = matrix(0, nrow=nm, ncol=nlr)


  for (j in 1:nlr) {
    cat(j, "of", nlr, "\n")
    lr = lib_ratios[j]
    for (i in 1:nm) {
      x = raw_expr_means[i]
      Xstats = compute_logX_distribution(x, dispersion, lr)
      logX_mean[i,j] = Xstats$logX_mean
      logX_var[i,j] = Xstats$logX_var
    }
  }

  saveRDS(list(logX_mean=logX_mean,
               logX_var=logX_var,
               lib_ratios=lib_ratios,
               raw_expr_means=raw_expr_means), rds_file)
  return (sdir)
}

compute_logX_distribution = function(xij, dispersion, library_ratio)
{
  if (xij==0)
    return (list(logX_mean=0, logX_var=0))

  max_val = 2*round(xij) + 1
  while (pnbinom(max_val, mu=xij, size=dispersion) < .99)
    max_val = 2*max_val
  vals = 0:max_val
  d = dnbinom(vals, mu=xij, size=dispersion)
  term1 = log1p(vals/library_ratio)
  term2 = d*term1 # d*log(1+vals/Lj*L)
  term3 = term2*term1  # d*log(1+vals/Lj*L)^2
  logX_mean = sum(term2)
  logX_var = sum(term3) - logX_mean^2

  return (list(logX_mean=logX_mean, logX_var=logX_var))
}

read_spike_processing = function(dispersion)
{
  sdir = file.path("spike_preprocessing")
  rds_file = file.path(sdir,
                       paste("spike_pre_dispersion_",  dispersion, ".rds",
                             sep=""))

  return (readRDS(rds_file))
}

