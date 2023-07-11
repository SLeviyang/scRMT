# Function below is useful for interactive analysis of the bulk.  Not
# used in the manuscript

#' Compare eigenvectors of expression matrix to marchenko-pastur
#' distribution.
#'
#' @param s Seurat object
#' @param X expression matrix, cells by genes
#' @param drop number of the dominant eigenvectors of X to drop.
#' Useful in examining bulk of X
#' @param scale_drop Should the remaining eigenvectors be rescaled
#' so their sum is the Frobenius norm of X?
#' @param binwidth bin size for histogram
#' @param permute should the gene expressions be permuted
#' independently for each gene to remove all correlations?
show_mp_approximation = function(s=NULL,
                                 X=NULL,
                                 drop=0,
                                 binwidth=.1,
                                 xlim=NULL,
                                 permute=F,
                                 include_limit=F)
{
  if (is.null(s) & is.null(X)) {
    message("Either s or X must be provided")
  }
  if (is.null(X))
    X = scaled_counts.dataset(s)

  if (drop > 0) {
    info = X_filter_eigenvalues(X, drop=drop, return_eigenvalues = T)
    X = info$X
    eX = info$eX
  }

  if (permute)
    X = apply(X, 2, function(g) g[sample.int(length(g))])

  if ((drop > 0 & permute) | drop == 0)
    eX = eigen(t(X) %*% X)$values

  df_values = data.frame(values=eX)
  # compute theoretical bulk
  if (include_limit) {
    p = param(colSums(X^2), w=1, gamma=ncol(X)/nrow(X))
    s = spectrode(p, dx=0.001)
    if (p$gamma < 1)
      s$density = s$density/p$gamma

    df_limit = data.frame(x=s$grid, y=s$density)
    max_val = max(c(df_values$values+1, df_limit$x+1))
  } else
    max_val = max(df_values$values+1)

  g = ggplot() +
    geom_histogram(data=df_values, aes(x=values),
                   binwidth=binwidth,
                   colour = 1, fill = "white")
  if (include_limit)
    g = g + geom_line(data=df_limit, aes(x=x,y=y*binwidth*length(eX)),
              color="red", size=1)
  if (is.null(xlim))
    g = g + xlim(0, max_val)
  else
    g = g + xlim(xlim)


  return (g)
}




# drops dominant eigenvectors
X_filter_eigenvalues = function(X, drop=NULL,
                                keep=NULL,
                                return_eigenvalues=F)
{
  if (is.null(keep) & is.null(drop)) {
    message("either drop or keep must not be NULL")
    return (NULL)
  }
  if (!is.null(drop)) {
    cat("dropping", drop, "eigenvalues", "\n")
    keep_vec = (drop+1):ncol(X)
  } else {
    cat("keeping", drop, "eigenvalues", "\n")
    keep_vec = 1:keep
  }
  
  e = eigen(t(X) %*% X)
  V = e$vectors[,keep_vec]
  X = X %*% V %*% t(V)
  eX = e$values[keep_vec]
  
  if (return_eigenvalues)
    return (list(X=X, eX=eX))
  else
    return (X)
}
