#' Make the filtering and clustering tables used in the manuscript
make.manuscript_DE = function(all_datasets=NULL)
{
  if (is.null(all_datasets))
    all_datasets = get_all_DE_datasets.dataset_DE()

  for (dataset in all_datasets) {
    if (dataset == "Intestine")
      permutation_types = c("none", "module")
    else
      permutation_types = "none"
    bulk_profile_table.dataset_DE(dataset, permutation_types)
    filtering_table.dataset_DE(dataset)
    filtering_theory_table.dataset_DE(dataset, permutation_types)
    clustering_table.manuscript_DE(dataset)
    clustering_table_theory.manuscript_DE(dataset, permutation_types)
  }
}

get_spike_k.manuscript_DE = function(dataset)
{
  switch(dataset,
         Zheng=40,
         HumanNasal=40,
         Intestine=40,
         stop("problem"))
}



###################

bulk_profile_table.dataset_DE = function(dataset,
                                         permutation_types="none",
                                         overwrite=F,
                                         max_spike=10,
                                         d_spike=0.001)
{
  t = here::here(table_dir.dataset_DE(dataset),
                 paste(dataset, "_bulk_profile_table.csv", sep=""))

  d = data.frame(permutation_type=permutation_types, stringsAsFactors = F)

  if (!file.exists(t) | overwrite) {
    d_list = parallel::mclapply(1:nrow(d), function(i) {

      pt = d[i,]
      cat(i, "of", nrow(d), "\n")
      # load the baseline
      s0 = nISG.dataset_DE(0, pt, dataset)

      # remove the leading singular values and compute the bulk
      uv = s0@misc$uv
      spike_k = switch(pt, single=3, get_spike_k.manuscript_DE(dataset))
      Xs = uv$u[,1:spike_k] %*% diag(uv$d[1:spike_k]) %*% t(uv$v[,1:spike_k])
      b = bulk.analysis_util(s0, Xs)
      
      # for different 
      spikes = seq(0, max_spike, d_spike)

      spike_p = spike_profile.analysis_util(b, spikes=spikes)
      spike = spike_p$spike
      sv = spike_p$singular_value
      u_proj2 = spike_p$u_proj2
      v_proj2 = spike_p$v_proj2

      d_out = data.frame(spike=spike,
                         singular_value=sv,
                         u_proj2=u_proj2,
                         v_proj2=v_proj2) %>%
        dplyr::mutate(dataset=dataset,
                      permutation_type=pt,
                      .before=1)
      return (d_out)
    }, mc.cores=3)
    d = do.call(rbind, d_list)
    write.csv(d, t, row.names = F)
  }

  d = read.csv(t, check.names=F)
  return (d)
}

######################################
# filtering

filtering_theory_table.dataset_DE = function(dataset,
                                             permutation_types="none",
                                                       overwrite=F,
                                                       A_max=10,
                                                       debug=F)
{
  t = here::here(table_dir.dataset_DE(dataset),
                 paste(dataset, "_filtering_theory.csv", sep=""))

  if (!file.exists(t) | overwrite) {

    spike = seq(0, sqrt(A_max), 0.01)
    d_spike = bulk_profile_table.dataset_DE(dataset)
    d = data.frame(permutation_type = permutation_types, stringsAsFactors = F)

    d_list = lapply(1:nrow(d), function(i) {

      pt = d[i,]
      cat(i, "of", nrow(d), "\n")


      if (debug) {
        b = filtering_table.dataset_DE(dataset) %>%
          dplyr::filter(permutation_type==pt, A_spike<=max(spike^2))
        spike = seq(0, sqrt(max(b$A)), .1)

      }

      s0 = nISG.dataset_DE(0, pt, dataset)


      # load the spike/bulk information,
      spike_p = dplyr::filter(d_spike, permutation_type==pt)
      spike_interpolation = bulk_profile_interpolation.analysis_util(spike_p$spike,
                                                                     spike_p$singular_value,
                                                                     spike_p$u_proj2)
      sv_f = spike_interpolation$singular_value
      coherent_proj2_f = spike_interpolation$coherent_proj2

      # compute the svd perturbation spike*u*v + \sum_i d_i x_i y_i^T
      rank = switch(pt, single=3, get_spike_k.manuscript_DE(dataset))
      Xs_uv = s0@misc$uv

      u = make_u_spike.analysis.util(s0$stim)
      sp = svd_perturbation.analysis_util(u,
                                          d=Xs_uv$d[1:rank],
                                          Xs_uv$u[,1:rank])

      if (pt == "single") {
        A_captured = (sv_f(spike)^2)*coherent_proj2_f(spike)
      } else if (pt == "none" | pt=="orthogonal" | pt=="module") {
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
        dplyr::mutate(dataset=dataset,
                      permutation_type=pt,
                      .before=1)

      if (debug) {

        plot(d_out$A, d_out$A_spike, type="l", ylim=c(0,max(c(d_out$A_spike,
                                                              b$A_spike))))
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

filtering_table.dataset_DE = function(dataset,
                                         overwrite=F,
                                                debug=F)
{
  ft_filename = here::here(table_dir.dataset_DE(dataset),
                           paste(dataset, "_filtering.csv", sep=""))

  d = get_generated_table.dataset_DE(dataset)

  if (nrow(d)==0) {
    cat("no data has not been generated\n")
    return (NULL)
  }

  if (debug) {
    d = dplyr::filter(d, permutation_type=="single")
  }


  if (!file.exists(ft_filename) | overwrite) {

    tb_full = parallel::mclapply(1:nrow(d), function(i) {
      cat(i, "of", nrow(d), "\n")

      nISG = d$nISG[i]
      pt = d$permutation_type[i]
      s = nISG.dataset_DE(nISG, pt, dataset)

      spike_k = bulk_k = switch(pt, single=3, get_spike_k.manuscript_DE(dataset))


      snr = s@misc$snr
      A = sum(snr/(1+snr))
      u_bar = make_u_spike.analysis.util(s$stim)
      df = filtering_measure.analysis_util(s, u_bar, spike_k, bulk_k)

      df = dplyr::mutate(df,
                         dataset=dataset,
                         permutation_type=pt,
                         nISG=nISG,
                         A=A, .before=1)

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


#####################################################
# clustering
#' Computes fraction of simulated cells that have a stimulated cell as
#' a nearest neighbor.  Theoretical fraction is also computed.
clustering_table.manuscript_DE = function(dataset,
                                          overwrite=F,
                                                  debug=F)
{
  ft_filename = here::here(table_dir.dataset_DE(dataset),
                 paste(dataset, "_clustering.csv", sep=""))

  d = get_generated_table.dataset_DE() %>%
    dplyr::filter(dataset==!!dataset)

  if (debug) {
    d = dplyr::filter(d, permutation_type=="single")
  }


  if (nrow(d)==0) {
    cat("data has not been generated\n")
    return (NULL)
  }

  if (!file.exists(ft_filename) | overwrite) {

    tb_full = lapply(1:nrow(d), function(i) {

      cat(i, "of", nrow(d), "\n")
      dd = d[i,]

      s = nISG.dataset_DE(dd$nISG, dd$permutation_type, dataset)
      snr = s@misc$snr
      A = sum(snr/(1+snr))

      pt = dd$permutation_type
      spike_k = switch(pt, single=3, get_spike_k.manuscript_DE(dataset))
      nn = nn_sampled.analysis_util(s, spike_k)

      df = d[i,] %>% dplyr::mutate(A=A,
                                   spike_k=spike_k,
                                   nn=nn)
      if (debug)
        browser()

      return(df)
    })

    tb = do.call(rbind, tb_full)
    write.csv(tb, ft_filename, row.names = F)
  }

  tb = read.csv(ft_filename, check.names = F)
  return (tb)
}


clustering_table_theory.manuscript_DE = function(dataset,
                                                 permutation_types="none",
                                                 overwrite=F,
                                                 debug=F)

{

  ft_filename = here::here(table_dir.dataset_DE(dataset),
                           paste(dataset, "_clustering_theory.csv", sep=""))

  d = data.frame(permutation_type=permutation_types,
                 stringsAsFactors = F)

  spike_p_all = bulk_profile_table.dataset_DE(dataset)
  sv_values = seq(0,sqrt(10), .05)


  if (!file.exists(ft_filename) | overwrite) {

    tb_full = lapply(1:nrow(d), function(i) {
      cat(i, "of", nrow(d), "\n")

      pt = d[i,]
      spike_k = bulk_k = switch(pt, single=3, get_spike_k.manuscript_DE(dataset))
      # debug
      if (debug) {
        b = clustering_table.manuscript_DE(dataset) %>%
            dplyr::filter(permutation_type==pt)
        sv_values = seq(0,sqrt(6), .05)
      }

      s0 = nISG.dataset_DE(0, pt, dataset)
      N = 20
     
      # bulk information
      spike_p = spike_p_all %>%
        dplyr::filter(permutation_type==pt)
      spike_interpolation = bulk_profile_interpolation.analysis_util(spike_p$spike,
                                                                     spike_p$singular_value,
                                                                     spike_p$u_proj2)

      # due to interpolation sometimes values slightly less than 0
      sv = pmax(0,spike_interpolation$singular_value(sv_values))
      alpha2 = pmax(0,spike_interpolation$coherent_proj2(sv_values))

      if (pt == "single") {
        nn = nn_theory_homegeneous.analysis_util(s0, spike_k, sv, alpha2,
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
        dplyr::mutate(dataset=dataset,
                      permutation_type=pt,
                      .before=1)

      if (debug) {
        plot(df$A, df$nn, type="l", ylim=c(0,1))
        points(b$A, b$nn, col="red")
        browser()
      }

      return(df)
    })

    tb = do.call(rbind, tb_full)
    write.csv(tb, ft_filename, row.names = F)
  }

  tb = read.csv(ft_filename, check.names = F)
  return (tb)
}


