#' Makes true and predicted clustering tables for the Kang et al dataset and
#' splatter dataset


make.manuscript_clustering = function(overwrite=F)
{
  clustering_table.manuscript_clustering(overwrite=overwrite)
  for (pt in c("none"))
    clustering_table_theory.manuscript_clustering(pt, overwrite=overwrite)
}


#' Computes theoretical fraction of simulated cells that have a stimulated cell as
#' a nearest neighbor. 
clustering_table_theory.manuscript_clustering = function(permutation_type,
                                                         overwrite=F,
                                                         debug=F)

{
  ft_filename = here::here(dir_analysis.Kang_IFN(),
                           paste("manuscript_clustering_theory_",
                                 permutation_type,
                                 ".csv", sep=""))


  d = get_generated_table.Kang_IFN() %>%
    dplyr::filter(!is.element(celltype,
                              c("CD4 T cells", "ct5"))) %>%
    dplyr::select(celltype, insilico, p, permutation_type, ISG_type) %>%
    dplyr::filter(permutation_type==!!permutation_type) %>%
    unique.data.frame()


  if (debug) {
    d = get_generated_table.Kang_IFN() %>%
        dplyr::select(celltype, insilico, p, permutation_type, ISG_type) %>%
        unique.data.frame()

    d = dplyr::filter(d, celltype=="Kang_all", p==0.5,
                      permutation_type=="none",
                    insilico==F)
  }

  spike_p_all = bulk_profile_table.analysis_util()
  sv_values = seq(0,sqrt(10), .05)


  if (nrow(d)==0) {
    cat("data has not been generated\n")
    return (NULL)
  }

  if (!file.exists(ft_filename) | overwrite) {

    tb_full = parallel::mclapply(1:nrow(d), function(i) {
      cat(i, "of", nrow(d), "\n")

      dd = d[i,]
      spike_k = bulk_k = get_spike_k.Kang_IFN(dd$celltype)


      # debug
      if (debug) {
        b = clustering_table.manuscript_clustering() %>%
          dplyr::filter(celltype==dd$celltype,
                        permutation_type==dd$permutation_type,
                        p==dd$p)
        sv_values = seq(0,sqrt(6), .05)
        #sv_values = seq(0,sqrt(max(b$A)), .05)
      }

      if (dd$celltype != "splatter_all" & dd$celltype != "Kang_all")
        ISG_type = NULL
      else
        ISG_type = d$ISG_type[i]

      s0 = nISG.Kang_IFN(0, dd$permutation_type,
                         dd$insilico, dd$celltype,
                         dd$p, ISG_type = ISG_type)
      # sample more for small cell types
      # if (ncol(s0) < 2500)
      #   N = 5
      # else
      #   N = 20
      N = 20

      # bulk information
      spike_p = spike_p_all %>%
        dplyr::filter(celltype==dd$celltype,
                      insilico==dd$insilico,
                      p==dd$p,
                      permutation_type==dd$permutation_type,
                      ISG_type==dd$ISG_type)
      spike_interpolation = bulk_profile_interpolation.analysis_util(spike_p$spike,
                                                                     spike_p$singular_value,
                                                                     spike_p$u_proj2)

      # due to interpolation sometimes values slightly less than 0
      sv = pmax(0,spike_interpolation$singular_value(sv_values))
      alpha2 = pmax(0,spike_interpolation$coherent_proj2(sv_values))

      if (dd$permutation_type == "single") {
        nn = nn_theory_homogeneous.analysis_util(s0, spike_k, sv, alpha2,
                                          N=N, debug=debug)
      }
      else {
        nn = nn_theory_none.analysis_util(s0@misc$uv$u,
                                          s0@misc$uv$d,
                                          s0$stim,
                                          sv_values, spike_k,
                                          spike_interpolation$singular_value,
                                          spike_interpolation$coherent_proj2,
                                          unlocalized=F, N=N,
                                          debug=debug)
      }


      df = data.frame(A=sv_values^2,
                      spike_A=sv^2,
                      alpha2=alpha2,
                      nn=nn$nn,
                      s2=nn$s2) %>%
            dplyr::mutate(celltype=dd$celltype,
                          insilico=dd$insilico,
                          p=dd$p,
                          permutation_type=dd$permutation_type,
                          ISG_type=dd$ISG_type, .before=1)

      if (debug) {
        plot(df$A, df$nn, type="l", ylim=c(0,1))
        #lines(df$A, df$nn_u, col="green")
        #lines(df$A, df$nn_s, col="blue")
        points(b$A, b$nn, col="red")
        browser()
      }

      return(df)
    }, mc.cores=3)

    tb = do.call(rbind, tb_full)
    write.csv(tb, ft_filename, row.names = F)
  }

  tb = read.csv(ft_filename, check.names = F)
  return (tb)
}


#' Computes fraction of simulated cells that have a stimulated cell as
#' a nearest neighbor.  Theoretical fraction is also computed.
clustering_table.manuscript_clustering = function(overwrite=F,
                                                debug=F)
{
  ft_filename = here::here(dir_analysis.Kang_IFN(),
                           "manuscript_clustering.csv")


  d = get_generated_table.Kang_IFN()

  if (nrow(d)==0) {
    cat("data has not been generated\n")
    return (NULL)
  }

  if (!file.exists(ft_filename) | overwrite) {

    tb_full = parallel::mclapply(1:nrow(d), function(i) {
      cat(i, "of", nrow(d), "\n")
      dd = d[i,]

      if (dd$celltype != "splatter_all" & dd$celltype != "Kang_all")
        ISG_type = NULL
      else
        ISG_type = d$ISG_type[i]

      s = nISG.Kang_IFN(dd$nISG, dd$permutation_type,
                        dd$insilico, dd$celltype,
                        dd$p, ISG_type)
      snr = s@misc$snr
      A = sum(snr/(1+snr))

      spike_k = get_spike_k.Kang_IFN(dd$celltype)
      nn = nn_sampled.analysis_util(s, spike_k)

      df = d[i,] %>% dplyr::mutate(A=A,
                                   spike_k=spike_k,
                                   nn=nn)
      if (debug)
        browser()

      return(df)
    }, mc.cores=3)

    tb = do.call(rbind, tb_full)
    write.csv(tb, ft_filename, row.names = F)
  }

  tb = read.csv(ft_filename, check.names = F)
  return (tb)
}
