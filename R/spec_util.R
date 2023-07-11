#' Plot the ESD of a spectrode object
#'
#' @param s a spectrode object
#' @param evals a numeric vector of eigenvalues that are to
#' be compared to the ESD.  Ignored if NULL
#' @param nbins number of bins for evals histogram.  Not relevant if
#' evals is NULL
#' @param remove_zero_evals should the 0's in evals be removed?   This
#' allows for appropriate comparison to the spectrode density, which does
#' not include 0's.  Only relevant when gamma < 1.
#' @param xlim a numeric vector of length 2 giving the endpoints of the
#' x-axis in the generated graph.  If NULL then defaults to including
#' all spectrum values.  Useful for zooming in on a particular portion
#' of the spectrum.
#'
#' @details If evals is NULL, the ESD is plotted.  Otherwise, evals
#' are shown in a histogram and ESD is plotted on top, adjusted for the
#' histogram binwidth.
#'
#' The spectrode density does not include zero eigenvalues so the plot defaults
#' to removing 0's in the passed evals to make the comparison appropriate
plot.spectrode = function(s,
                          evals=NULL,
                          nbins=30,
                          xlim=NULL,
                          binwidth=NULL,
                          remove_zero_evals=T)
{
  param = s$param
  df_theory = data.frame(x=s$grid,
                         y=s$density)
  if (param$gamma < 1)
    df_theory$y = df_theory$y/param$gamma

  if (is.null(evals))
    g = ggplot() +
    geom_line(data=df_theory, aes(x=x,y=y),
              color="red", size=1)
  else {
    if (is.null(binwidth))
      binwidth = (max(evals) - min(evals))/nbins
    df_true = data.frame(x=evals)
    g = ggplot() +
      geom_histogram(data=df_true, aes(x=x),
                     binwidth=binwidth,
                     colour = 1, fill = "white") +
      geom_line(data=df_theory, aes(x=x,y=y*binwidth*length(evals)),
                color="red", size=1)

  }


  if (is.null(xlim))
    g = g + xlim(0, max(evals, s$grid))
  else
    g = g + xlim(xlim)

  return (g)
}


#' Plots the spectrum of a finite sample matrix and the limiting
#' spectrum as computed by spectrode.
#'
#' Useful in assessing the accuracy of spectrode to finite samples under
#' the model assumptions
#'
#' @param param sprectrodeR parameter object
#' @param n number of rows in the matrix.  The param object
#' then specifies the number of cols based on the gamma value.
#' @param dx grid size at which the limiting spectrum density is computed
#' @param nbins, number of bins to use in showing the finite sample
#' spectrum distribution
#'
#' @returns a ggplot object
#' @export
simulate.spectrode = function(param, n=2000,
                             dx=0.01, nbins=30,
                             xlim=NULL)
{
  t = param$t; w = param$w; gamma = param$gamma
  nt = length(t)
  if (length(w) != nt) {
    message("invalid parameters.  w must have same length as t")
    return (NULL)
  }

  if (abs(sum(w) - 1) > 1E-8) {
    message("w must be normalized to sum to 1")
    return (NULL)
  }

  # simulate
  N = round(n*gamma)
  Ni = round(w*N)

  Xl = lapply(1:nt, function(i) {
    matrix(rnorm(n*Ni[i]), nrow=n, ncol=Ni[i])/sqrt(n)*sqrt(t[i])
  })
  X = do.call(cbind, Xl)

  if (gamma < 1)
    e = eigen(t(X) %*% X)$values
  else
    e = eigen(X %*% t(X))$values

  # theory
  s = spectrode(param, dx=dx)
  g = plot.spectrode(s, evals=e,
                     nbins=nbins,
                     xlim=xlim,
                     remove_zero_evals=T)

  return (g)
}

