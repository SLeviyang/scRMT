# Figures that do not fit into the bulk/filtering/clustering topics

make.manuscript_odds_and_ends = function(m_version)
{
  fname = here::here(figure_dir.manuscript(m_version),
                     "freq_filter_cluster.pdf")
  if (!file.exists(fname)) {
    g = frequency_effects.manuscript_odds_ends(max_A=4)
    pdf(fname, width=6, height=4)
    print(g)
    dev.off()
  }


  fname = here::here(figure_dir.manuscript(m_version),
                     "BBP_example.pdf")
  if (!file.exists(fname)) {
    pdf(fname, width=8, height=1.5)
    g = BBP_example.manuscript_odds_and_ends() +
    theme(strip.text.x = element_text(size = 12, face="bold")) +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=14,face="bold"))
    print(g)
    dev.off()
  }

}

#' This figure shows effects of changing p_1 on filtering and clustering
frequency_effects.manuscript_odds_ends = function(max_A=4,
                                                  celltype="NK cells",
                                                  pt="none")
{
  p_vals = c(0.05, 0.25, 0.50)

  # nearest neighbor plots
  df = clustering_table.manuscript_clustering() %>%
    dplyr::filter(celltype==!!celltype,
                  A <= max_A,
                  is.element(p, p_vals),
                  permutation_type==pt)
  dt = clustering_table_theory.manuscript_clustering(pt) %>%
    dplyr::filter(celltype==!!celltype,
                  A <= max_A,
                  is.element(p, p_vals))

  g2 = ggplot() + geom_point(aes(x=A, y=nn), data=df) +
    geom_line(aes(x=A, y=nn), data=dt) +
    facet_wrap("p", nrow=1)  +
    xlab("signal strength") + ylab("nn metric") +
    theme(strip.text.x = element_text(size = 12, face="bold")) +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=14,face="bold")) +
    scale_y_continuous(breaks = c(0,.5,1),
                       labels=c("0", ".5", "1"))

  return (g2)

}



#' An example of a weak ISG referred to in the text, no figure
weak_ISG_example.manuscript_odds_ends = function()
{
  s = permutation.Kang_IFN("none", F, "CD14+ Monocytes", 0.5)

  Y = scaled_counts.dataset(s)
  D = s@assays$RNA@data %>% Matrix::as.matrix()

  u = make_u_spike.analysis.util(s$stim)
  cc = ((t(Y) %*% u)^2)
  delta_mu = 2*(D %*% u)/sqrt(nrow(Y))
  ISG = ifelse(is.element(colnames(Y), s@misc$ISG), "ISG", "baseline")

  df = data.frame(cc=cc, gene=colnames(Y),
                  delta_mu=delta_mu,
                    ISG=ISG,
                    celltype=s$celltype[1]) %>%
      dplyr::arrange(desc(cc))

  ind_nonISG = which(df$ISG != "ISG")[1]
  print(df[ind_nonISG,])
  cat("median correlation of ISG",
      dplyr::filter(df, ISG=="ISG") %>%
        pull(cc) %>% median, "\n")
  print(head(df))

  return (NULL)
}
  


#' This figure shows the hat(lambda) and alpha of the BBP critical threshold
BBP_example.manuscript_odds_and_ends = function()
{
  d_spike = bulk_profile_table.analysis_util()

  # load the spike/bulk information,
  spike_p = dplyr::filter(d_spike, celltype=="Kang_all",
                          insilico==F,
                          p==0.5,
                          permutation_type=="none")
  spike_interpolation = bulk_profile_interpolation.analysis_util(spike_p$spike,
                                                                 spike_p$singular_value,
                                                                 spike_p$u_proj2)
  sv_f = spike_interpolation$singular_value
  coherent_proj2_f = spike_interpolation$coherent_proj2


  lambda = seq(0, 4, .01)
  lambda_hat = sv_f(lambda)
  # fix lambda_hat so that it's b when below BBP
  b = min(lambda_hat[lambda_hat > 1E-2])
  lambda_hat[lambda_hat < 1E-2] = b
  alpha2 = coherent_proj2_f(lambda)
  alpha2 = ifelse(alpha2 < 1E-10, 0, alpha2)
  tau = lambda[min(which(alpha2 > 1E-10))]
  cat("critical point", tau, "\n")

  df1 = data.frame(lambda=lambda, value=lambda_hat^2)
  g1 = ggplot(data=df1) + geom_line(aes(x=lambda, y=value)) +
       xlab("signal strength") +
       ylab(latex2exp::TeX("$\\hat{lambda}^2$")) +
       ylim(0, max(df1$value))

  df2 = data.frame(lambda=lambda, value=alpha2)
  g2 = ggplot(data=df2) + geom_line(aes(x=lambda, y=value)) +
       geom_vline(aes(xintercept=tau), color="red") +
       xlab("signal strength") +
       ylab(latex2exp::TeX("$\\alpha^2$")) + ylim(0,1)

  df3 = data.frame(lambda=lambda, value=lambda_hat*sqrt(alpha2))
  g3 = ggplot(data=df3) + geom_line(aes(x=lambda, y=value)) +
   ylim(0, max(df1$value)) +
    xlab("signal strength") +
    ylab("filtered ss")

  g = ggpubr::ggarrange(plotlist=list(g1, g2),
                    ncol=2, nrow=1)
  return (g)
}

#' Correlations results referred to in the text
celltype_baseline_correlations = function()
{
  cc = plyr::adply(get_celltypes.Kang(), 1,  function(celltype) {

    if (celltype != "Kang_all")
      s = nISG.Kang_IFN(0, "none", F, celltype, p=0.5)
    else
      s = nISG.Kang_IFN(0, "none", F, celltype, p=0.5, "joint")
    x = s@misc$uv$u
    u = make_u_spike.analysis.util(s$stim)

    w = replicate(20, {
      xx = apply(x, 2, function(z) z[sample.int(length(z))])
      max(abs(t(xx) %*% u))
    }) %>% max



    df = data.frame(dataset=celltype,
                    corr=max(abs(t(x) %*% u)),
                    rand=w)

    return (df)
  })

  return (cc)
}

