make.manuscript_filtering = function(overwrite=F)
{
  filtering_table.manuscript_filtering(overwrite=overwrite)
  filtering_theory_table.manuscript_filtering(overwrite=overwrite)
}


##################################################################
# data generation functions

#' Computes spike strength assuming bulk is a Marchenko-Pastur distribution
#'
#' @details Use the Kang_IFN Seurat object with n_ISG==0 to compute a bulk.
#' Based on that bulk, a spectrode object is created to compute the ESD
#' and a spike profile is produced.
filtering_theory_table.manuscript_filtering = function(overwrite=F,
                                                       A_max=6,
                                                       debug=F)
{
  t = here::here(dir_analysis.Kang_IFN(),
            "manuscript_filtering_theory.csv")

  if (!file.exists(t) | overwrite) {

    d = get_generated_table.Kang_IFN() %>%
        dplyr::select(celltype, insilico, p, permutation_type, ISG_type) %>%
        unique.data.frame()


    spike = seq(0, sqrt(A_max), 0.01)
    d_spike = bulk_profile_table.analysis_util(overwrite=F)


    if (debug) {
      d = dplyr::filter(d, celltype=="Kang_all", insilico==F, p==0.5,
                        permutation_type=="single")
    }

    d_list = lapply(1:nrow(d), function(i) {

      dd = d[i,]
      cat(i, "of", nrow(d), "\n")

      if (dd$ISG_type == "cellspecific")
        It = NULL
      else
        It = dd$ISG_type

      if (debug) {
        b = filtering_table.manuscript_filtering() %>%
          dplyr::filter(celltype==dd$celltype,
                        permutation_type==dd$permutation_type,
                        p==dd$p, A_spike<=max(spike^2))
        spike = seq(0, sqrt(max(b$A)), .1)

      }

      s0 = nISG.Kang_IFN(0, dd$permutation_type, dd$insilico,
                         dd$celltype, dd$p, It)

      # load the spike/bulk information,
      spike_p = dplyr::filter(d_spike, celltype==dd$celltype,
                                insilico==dd$insilico,
                                p==dd$p,
                                permutation_type==dd$permutation_type,
                                ISG_type==dd$ISG_type)
      spike_interpolation = bulk_profile_interpolation.analysis_util(spike_p$spike,
                                                                     spike_p$singular_value,
                                                                     spike_p$u_proj2)
      sv_f = spike_interpolation$singular_value
      coherent_proj2_f = spike_interpolation$coherent_proj2

      # compute the svd perturbation spike*u*v + \sum_i d_i x_i y_i^T
      rank = get_spike_k.Kang_IFN(dd$celltype)
      Xs_uv = s0@misc$uv

      u = make_u_spike.analysis.util(s0$stim)
      sp = svd_perturbation.analysis_util(u,
                                          d=Xs_uv$d[1:rank],
                                          Xs_uv$u[,1:rank])

      pt = dd$permutation_type
      if (pt == "single") {
        A_captured = (sv_f(spike)^2)*coherent_proj2_f(spike)
      } else if (pt == "none" | pt=="orthogonal") {
        A_captured = sapply(spike, function(c_spike) {

          # get the perturbed svd
          svd_p = compute_svd_perturbation.analysis_util(c_spike, sp)
          rho = sqrt(svd_p$eigenvalues)
          beta2 = svd_p$u_proj2

          # run the singular values through the bulk
          tilde_rho = sv_f(rho)
          # portion of the singular vectors that are not noise/coherent
          alpha2 = coherent_proj2_f(rho)
          captured = sum(alpha2*tilde_rho^2*beta2)


          return (captured)
        })
      } else {
        stop("permutation type not recognized!")
      }


      d_out = data.frame(A=spike^2,A_spike=A_captured) %>%
              dplyr::mutate(celltype=dd$celltype,
                            insilico=dd$insilico,
                            p=dd$p,
                            permutation_type=dd$permutation_type,
                            ISG_type=dd$ISG_type,
                            .before=1)

      if (debug) {

        plot(d_out$A, d_out$A_spike, type="l", ylim=c(0,max(c(d_out$A_spike,
                                                              b$A_spike))))
        #lines(df$A, df$nn_u, col="green")
        #lines(df$A, df$nn_s, col="blue")
        points(b$A, b$A_spike, col="red")
        browser()
      }

      return (d_out)
    })
    d = do.call(rbind, d_list)
    write.csv(d, t, row.names = F)
  }

  d = read.csv(t, check.names=F)
  return (d)
}

filtering_table.manuscript_filtering = function(overwrite=F,
                                                debug=F)
{
  ft_filename = here::here(dir_analysis.Kang_IFN(),
                           "manuscript_filtering.csv")

  d = get_generated_table.Kang_IFN()

  if (nrow(d)==0) {
    cat("no data has not been generated\n")
    return (NULL)
  }

  if (!file.exists(ft_filename) | overwrite) {

    tb_full = parallel::mclapply(1:nrow(d), function(i) {
      cat(i, "of", nrow(d), "\n")

      nISG = d$nISG[i]
      celltype = d$celltype[i]
      insilico = d$insilico[i]

      if (celltype != "splatter_all" & celltype != "Kang_all")
        ISG_type = NULL
      else
        ISG_type = d$ISG_type[i]
      s = nISG.Kang_IFN(nISG, d$permutation_type[i],
                        insilico, celltype,
                        d$p[i], ISG_type)

      spike_k = bulk_k = get_spike_k.Kang_IFN(celltype)

      snr = s@misc$snr
      A = sum(snr/(1+snr))
      u_bar = make_u_spike.analysis.util(s$stim)
      df = filtering_measure.analysis_util(s, u_bar, spike_k, bulk_k)

      df = dplyr::mutate(df,
                         celltype=celltype,
                         insilico=insilico,
                         p=d$p[i],
                         permutation_type=d$permutation_type[i],
                         ISG_type=d$ISG_type[i],
                         nISG=nISG,
                         A=A, .before=1)

      if (debug)
        browser()

      return(df)
    }, mc.cores=5)

    tb = do.call(rbind, tb_full)
    write.csv(tb, ft_filename, row.names = F)
  }

  tb = read.csv(ft_filename, check.names = F)
  return (tb)
}
