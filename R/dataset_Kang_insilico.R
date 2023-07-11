dir.Kang_insilico = function()
{
  d = here::here("data", "Kang_insilico_datasets")
  if (!dir.exists(d))
    dir.create(d)

  return (d)
}

get_celltypes.Kang_insilico = function(full=F)
{
  d = data.frame(n=c(500, 1000, 2000, 5000, 10000))
  d$celltype = paste("ct", 1:nrow(d), sep="")

  d = rbind(d, data.frame(celltype="splatter_all", n=sum(d$n)))
  if (full)
    return (d)
  else
    return (d$celltype)
}



ISG_table.Kang_insilico = function(celltype, overwrite=F)
{
  t_file = here::here(dir.Kang_insilico(),
                      paste("insilico_ISG_table_", celltype, ".csv", sep=""))

  if (!file.exists(t_file) | overwrite) {
    s = celltype.Kang_insilico(celltype)
    if (celltype != "splatter_all") {
      ISG_pattern = s@misc$patterns[paste("ISG_cellspecific_", celltype, sep="")]
      ISGs = which(s@misc$gene_patterns==ISG_pattern) %>% names
      d = data.frame(type="cellspecific", celltype=celltype,
                     gene=ISGs, stringsAsFactors = F)
    } else {

      gp = s@misc$gene_patterns
      ISG_joint_pattern = s@misc$patterns["ISG_joint"]
      ct = s$celltype %>% unique %>% sort

      d_joint = data.frame(type="joint", celltype=NA,
                           gene=which(gp==ISG_joint_pattern) %>% names,
                           stringsAsFactors = F)
      d_cell = plyr::adply(ct, 1, function(cct) {
        pat = s@misc$patterns[paste("ISG_cellspecific_", cct, sep="")]
        data.frame(type="cellspecific", celltype=cct,
                   gene=which(gp==pat) %>% names)
      }, .id = NULL)

      d = rbind(d_joint, d_cell)
    }

    write.csv(d, t_file, row.names = F)
  }

  d = read.csv(t_file, header=T, check.names = F)
  return (d)
}

make.Kang_insilico = function()
{
  make_celltype.Kang_insilico()
  make_all.Kang_insilico()
}

make_all.Kang_insilico = function()
{
   celltype_all.Kang_insilico()
   processed_celltype.Kang_insilico("splatter_all")
   stim_freq.Kang_insilico("splatter_all", 0.5)
   ISG_table.Kang_insilico("splatter_all")
}

#' Make all the ct, but not all
make_celltype.Kang_insilico = function()
{
  ct = get_celltypes.Kang_insilico() %>% setdiff("splatter_all")
  parallel::mclapply(ct, function(celltype) {
    cat("making cell type file", celltype, "\n")
    rds_file = here::here(dir.Kang_insilico(),
                          paste("insilico_", celltype, ".rds",sep=""))


    if (!file.exists(rds_file))
      celltype.Kang_insilico(celltype)
    return (NULL)
  }, mc.cores=5)

  parallel::mclapply(ct, function(celltype) {
    cat("making processed cell type file", celltype, "\n")
    rds_file = here::here(dir.Kang_insilico(),
                          paste("processed_insilico_", celltype, ".rds",sep=""))


    if (!file.exists(rds_file))
      processed_celltype.Kang_insilico(celltype)
    return (NULL)
  }, mc.cores=5)


  p = c(0.5)
  d = expand.grid(celltype=ct, p=p, stringsAsFactors = F)

  parallel::mclapply(1:nrow(d), function(i) {
    print(d[i,])

    celltype = d$celltype[i]
    p = d$p[i]
    rds_file = here::here(dir.Kang_insilico(),
                          paste("insilico_", celltype, "_", p, ".rds",sep=""))

    if (!file.exists(rds_file))
      s = stim_freq.Kang_insilico(celltype, p)

    return (NULL)
  }, mc.cores=5)

  for (cct in ct)
    ISG_table.Kang_insilico(cct)
}

#################################################################
#' Filters simulated cells to make the stimulated
#' cells a frequency p of the total number of cells
#'
#' @param celltype
#' @param p num_stim/total
stim_freq.Kang_insilico = function(celltype, p, overwrite=F)
{
  if (p > 0.5) {
    message("p must not exceed 0.5")
    return (NULL)
  }

  rds_file = here::here(dir.Kang_insilico(),
                        paste("insilico_", celltype, "_", p, ".rds",sep=""))

  if (!file.exists(rds_file) | overwrite) {
    s = celltype.Kang_insilico(celltype)
    all_ct = s$celltype %>% unique %>% sort

    ind_list = lapply(all_ct, function(ct) {
      ind_stim = which(s$stim == "stim" & s$celltype == ct)
      ind_ctrl = which(s$stim == "ctrl" & s$celltype == ct)

      if (length(ind_ctrl) != length(ind_stim))
        stop("this shouldn't happen!")

      n = length(ind_stim)
      nn = round(n*p/(1-p))

      if (nn > 0)
        return (list(stim=ind_stim[1:nn], ctrl=ind_ctrl))
      else
        return (list(stim=NULL, ctrl=ind_ctrl))
    })
    ind_stim = lapply(ind_list, "[[", "stim") %>% unlist
    ind_ctrl = lapply(ind_list, "[[", "ctrl") %>% unlist

    s = s[,c(ind_ctrl, ind_stim)]
    s = update_Seurat_splatter_info(s) %>%
        recompute_spike_from_Seurat()

    saveRDS(s, rds_file)
  }

  a = readRDS(rds_file)
  return (a)
}

#' Performs Seurat workflow on the celltype.Kang output.  Also
#' computes an svd for downstream analysis.
processed_celltype.Kang_insilico = function(celltype, overwrite=F)
{
  rds_file = here::here(dir.Kang_insilico(),
                        paste("processed_insilico_", celltype, ".rds",sep=""))

  if (!file.exists(rds_file) | overwrite) {
    s = celltype.Kang_insilico(celltype) %>%
      process.dataset(nfeatures=1000, scale.max=Inf, npcs=20)

    s = update_Seurat_splatter_info(s) %>%
        recompute_spike_from_Seurat()

    Y = scaled_counts.dataset(s)
    s@misc$uv = irlba::svdr(Y, k=20)

    saveRDS(s, rds_file)
  }

  a = readRDS(rds_file)
  return (a)
}

# Dataset to match each of the Kang celltypes
celltype.Kang_insilico = function(celltype, overwrite=F)
{
  if (!is.element(celltype, get_celltypes.Kang_insilico())) {
    cat(celltype, "is not a valid cell type")
    return (NULL)
  }

  if (celltype == "splatter_all") {
    a = celltype_all.Kang_insilico(overwrite=overwrite)
    return (a)
  }
  rds_f = file.path(dir.Kang_insilico(),
                    paste("insilico_", celltype, ".rds", sep=""))

  if (!file.exists(rds_f) | overwrite) {
    ct = celltype
    n = get_celltypes.Kang_insilico(full=T) %>% dplyr::filter(celltype==ct) %>%
        dplyr::pull(n)
    stim_n = n/2
    ctrl_n = n/2
    ctrl_n_1 = round(2/3*ctrl_n)
    ctrl_n_2 = ctrl_n - ctrl_n_1
    ct_sizes = c(ctrl_n_1, ctrl_n_2, stim_n) %>%
               setNames(paste(celltype, "_", 1:3, sep=""))

    stim = c(rep("ctrl", ctrl_n), rep("stim", stim_n))
    # non_DE, shared DE, ISGs
    patterns = c("000", "011", "001") %>%
               setNames(c("homogeneous", "ISG_weak",
                          paste("ISG_cellspecific_", celltype, sep="")))
    nISGs = 100
    # I make 6 times the ISGs so that when p=0.05 I can get the
    # same spike strength as when p=.5
    ngenes = c(850, 50, 6*nISGs)

    s = splatter(ct_sizes, patterns, ngenes, dispersion=10,
                 de.mean=12/nISGs, de.scale=sqrt(5/nISGs))

    s@meta.data$stim = stim
    s@misc$patterns = patterns

    s@misc$celltype = NULL
    s@meta.data$true_celltype = s$celltype
    s@meta.data$celltype = rep(celltype, ncol(s))

    saveRDS(s, rds_f)
  }

  a = readRDS(rds_f)
  return (a)
}

celltype_all.Kang_insilico = function(overwrite=F)
{
  rds_f = file.path(dir.Kang_insilico(), "splatter_all.rds")

  if (!file.exists(rds_f) | overwrite) {
    ct_info = get_celltypes.Kang_insilico(full=T) %>%
              dplyr::filter(celltype != "splatter_all", celltype != "ct5")

    n_u_ct = ct_info$celltype %>% length

    # each celltype has splits into ctrl and stim
    ct_par = lapply(1:nrow(ct_info), function(i) {
      n = ct_info$n[i]
      ct = ct_info$celltype[i]

      stim_n = n/2
      ctrl_n = n/2
      sizes = c(ctrl_n,  stim_n)
      stim = c(rep("ctrl", ctrl_n), rep("stim", stim_n))
      cell_label = paste("ct", i, c("ctrl", "stim"), sep="")
      celltype = rep(paste("ct", i, sep=""), n)

      return (list(sizes=sizes, stim=stim,
                   cell_label=cell_label,
                   celltype=celltype))
    })

    ct_sizes = lapply(ct_par, "[[", "sizes") %>% unlist
    names(ct_sizes) = lapply(ct_par, "[[", "cell_label") %>% unlist
    stim =  lapply(ct_par, "[[", "stim") %>% unlist
    celltype = lapply(ct_par, "[[", "celltype") %>% unlist

    patterns = c(paste(rep("00", n_u_ct), collapse=""),
                 paste(rep("01", n_u_ct), collapse="")) %>%
              setNames(c("homogeneous", "ISG_joint"))
    ct_ISG_patterns = sapply(1:n_u_ct, function(i) {
      u_pre = paste(rep("0", 2*(i-1)), collapse="")
      u_mid = "01"
      u_post = paste(rep("0", 2*(n_u_ct-i)), collapse="")
      patt = paste(u_pre, u_mid, u_post, sep="")
    }) %>% setNames(paste("ISG_cellspecific_", ct_info$celltype, sep=""))

    ct_identity_patterns = sapply(1:n_u_ct, function(i) {
      u_pre = paste(rep("0", 2*(i-1)), collapse="")
      u_mid = "11"
      u_post = paste(rep("0", 2*(n_u_ct-i)), collapse="")
      patt = paste(u_pre, u_mid, u_post, sep="")
    }) %>% setNames(paste("cell_identity_", ct_info$celltype, sep=""))
    patterns = c(patterns, ct_ISG_patterns, ct_identity_patterns)

    nISGs = 200
    ngenes = c(homogeneous=300,
               joint=nISGs,
               cellspecific_ISG=rep(40, n_u_ct),
               cellspecific_ISG=rep(100, n_u_ct))

    s = splatter(ct_sizes, patterns, ngenes, dispersion=10,
                 de.mean=20/nISGs, de.scale=sqrt(5/nISGs))

    s@meta.data$stim = stim
    s@meta.data$celltype = celltype
    s@misc$patterns=patterns
    s@misc$celltype = NULL

    saveRDS(s, rds_f)
  }

  a = readRDS(rds_f)
  return (a)

}
