# functions to save/access bulk information for Kang et al datasets

# For each nISG compute the bulk profile
bulk_profile_table.analysis_util = function(overwrite=F,
                                            max_spike=10,
                                            d_spike=0.01)
{
  t = here::here(dir_analysis.Kang_IFN(),
                 "bulk_profile_table.csv")
  
  if (!file.exists(t) | overwrite) {
    
    d = get_generated_table.Kang_IFN() %>%
      dplyr::select(celltype, insilico, p, permutation_type, ISG_type) %>%
      unique.data.frame()
    
    d_list = parallel::mclapply(1:nrow(d), function(i) {
      
      dd = d[i,]
      cat(i, "of", nrow(d), "\n")
      if (dd$ISG_type == "cellspecific")
        It = NULL
      else
        It = dd$ISG_type
      s0 = nISG.Kang_IFN(0, dd$permutation_type, dd$insilico,
                         dd$celltype, dd$p, It)
      
      # in the case of real data, estimate the mean
      if (dd$insilico) {
        Xs = s0@misc$Xs %>% t
      } else {
        uv = s0@misc$uv
        Xs = uv$u %*% diag(uv$d) %*% t(uv$v)
      }
      
      b = bulk.analysis_util(s0, Xs)
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
        dplyr::mutate(celltype=dd$celltype,
                      insilico=dd$insilico,
                      p=dd$p,
                      permutation_type=dd$permutation_type,
                      ISG_type=dd$ISG_type,
                      .before=1)
      return (d_out)
    }, mc.cores=5)
    d = do.call(rbind, d_list)
    write.csv(d, t, row.names = F)
  }
  
  d = read.csv(t, check.names=F)
  return (d)
}

#' Extract singular value and alpha2 information from the bulk profile table
spike_from_bulk_profile_table.analysis_util = function(nISG,
                                                       permutation_type,
                                                       insilico,
                                                       celltype,
                                                       p,
                                                       ISG_type="cellspecific")
{
  # loading data
  s = nISG.Kang_IFN(nISG=nISG, permutation_type=permutation_type,
                    insilico=insilico, celltype=celltype, p=p)
  snr = s@misc$snr
  A = sum(snr/(1+snr))
  true_sv = sqrt(A)
  
  # loading spike information
  spike_p = bulk_profile_table.analysis_util() %>%
    dplyr::filter(celltype==!!celltype,
                  insilico==!!insilico,
                  p==!!p,
                  permutation_type==!!permutation_type,
                  ISG_type==!!ISG_type)
  spike_interpolation = bulk_profile_interpolation.analysis_util(spike_p$spike,
                                                                 spike_p$singular_value,
                                                                 spike_p$u_proj2)
  sv = max(0,spike_interpolation$singular_value(true_sv))
  alpha2 = max(0,spike_interpolation$coherent_proj2(true_sv))
  
  return (list(true_sv=true_sv, sv=sv, alpha2=alpha2))
}


