# functions to define dataset I/O

dir.dataset_DE = function(dataset)
{
  d = here::here("analysis", paste(dataset, "_DE_datasets", sep=""))
  if (!dir.exists(d))
    dir.create(d)

  return (d)
}

table_dir.dataset_DE = function(dataset)
{
  {
    d = here::here("analysis", paste(dataset, "_DE_tables", sep=""))
    if (!dir.exists(d))
      dir.create(d)

    return (d)
  }
}


get_unprocessed.dataset_DE = function(dataset)
{
  s = switch(dataset,
             Zheng = unprocessed.Zheng(),
             HumanNasal = unprocessed.HumanNasal(),
             Intestine = unprocessed.Intestine(),
             stop("dataset not recognized"))

  return (s)
}


get_all_DE_datasets.dataset_DE = function()
{
  return (c("Intestine", "HumanNasal", "Zheng"))
}

#####################################

get_generated_table.dataset_DE = function(all_datasets=NULL)
{
  if (is.null(all_datasets))
    all_datasets = get_all_DE_datasets.dataset_DE()

  tb = plyr::adply(all_datasets, 1, function(dataset) {
    files = dir(dir.dataset_DE(dataset))
    files = files[!grepl("permutation", files) & !grepl("base", files) & !grepl("gene_table", files)]
    files_no_rds = sapply(files, function(s) substr(s, 1, nchar(s)-4))

    ss = strsplit(files_no_rds, split="_")
    nISG = sapply(ss, "[", 2) %>% as.numeric
    pt = sapply(ss, "[", 3)

    df = data.frame(dataset=dataset,
                    nISG=nISG, permutation_type=pt, stringsAsFactors = F) %>%
      dplyr::arrange(nISG)

    return(df)
  })

  return (tb)
}




make.dataset_DE = function(all_datasets=NULL, N=40)
{
  if (is.null(all_datasets))
    all_datasets = get_all_DE_datasets.dataset_DE()

  for (dataset in all_datasets) {
    s = base.dataset_DE(dataset)
    parallel::mclapply(c("none"), function(pt) {
      sp = permutation.dataset_DE(pt, dataset)
      return (NULL)
    }, mc.cores = 3)

    snr = s@misc$snr
    sig = cumsum(snr/(1+snr))

    # get up to signal = 10
    max_BBP = max(which(sig <= 2))
    max_ind = max(which(sig <= 10))
    nISG_BBP = c(0,1,seq(0, max_BBP, length.out=N/2))
    nISG_big = seq(max_BBP, max_ind, length.out=N/2)

    nISG = c(nISG_BBP, nISG_big) %>% round %>% unique %>% sort

    parallel::mclapply(nISG, function(cISG) {
    #  nISG.dataset_DE(cISG, permutation_type="single", dataset)
    #  nISG.dataset_DE(cISG, permutation_type="orthogonal", dataset)
      nISG.dataset_DE(cISG, permutation_type="none", dataset)
      if (dataset=="Intestine")
        nISG.dataset_DE(cISG, permutation_type="module", dataset)
      return (NULL)
    }, mc.cores=3)
  }
}



##############################################################
#'
base.dataset_DE = function(dataset, overwrite=F)
{
  rds_file = file.path(dir.dataset_DE(dataset), paste(dataset, "_base.rds", sep=""))
  if (!file.exists(rds_file) | overwrite) {
    s = get_unprocessed.dataset_DE(dataset)
    gene_table = gene_table.dataset_DE(dataset)

    all_DE = intersect(rownames(s), gene_table$gene)
    non_DE = setdiff(rownames(s), all_DE)

    s = s[c(non_DE, all_DE),]
    s = process.dataset(s, nfeatures=1000, scale.max=Inf, npcs=20)

    active_DE = intersect(rownames(s), all_DE)
    s@misc$DE = active_DE

    # pick the snr that are greater than 0.1 as the head
    snr = snr.analysis_util(s, s$stim, s@misc$DE)
    ind = order(snr, decreasing=T)
    snr = snr[ind]
    active_cutoff = max(which(snr > .1))

    # to resolve the BBP make the initial gene snr 0.1
    if (active_cutoff < length(snr))
      shuffle_ind = c(rev(1:active_cutoff), (active_cutoff+1):length(snr))
    else
      shuffle_ind = rev(1:length(snr))
    snr = snr[shuffle_ind]
    s@misc$DE = names(snr)
    s@misc$snr = snr

    Y = scaled_counts.dataset(s)
    uv = irlba::svdr(Y[,s@misc$DE], k=40)
    rownames(uv$v) = s@misc$DE
    s@misc$uv = uv

    saveRDS(s, rds_file)
  }

  s = readRDS(rds_file)
  return (s)
}


permutation.dataset_DE = function(permutation_type, dataset, k=20,
                             overwrite=F)
{
  if (!is.element(permutation_type, c("single", "orthogonal", "none", "module"))) {
    message("invalid permutation type")
    return (NULL)
  }

  rds_file = here::here(dir.dataset_DE(dataset),
                        paste(dataset, "_", "permutation_", permutation_type, ".rds",
                              sep=""))

  if (!file.exists(rds_file) | overwrite) {
    s = base.dataset_DE(dataset)
    DE = is.element(rownames(s), s@misc$DE)
    ncells = ncol(s)

    # I always permute with classes in s$stim which means
    # that DE genes are kept B specific, but all other celltype
    # info is lost.
    per = switch(permutation_type,
                 module=permutation_module.dataset_util(DE, s$stim),
                 single=permutation_single.dataset_util(DE, s$stim),
                 orthogonal=permutation_orthogonal.dataset_util(ncells, DE),
                 none=permutation_none.dataset_util(ncells, DE),
                 stop("base permutation type"))

    Y = s@assays$RNA@counts %>% Matrix::as.matrix() %>%
      apply_permutation.Kang_IFN(per)

    stim_perm = per[which(DE)[1],]
    celltype_perm = per[which(!DE)[1],]

    s2 = Seurat::CreateSeuratObject(Y,
                                    project = dataset,
                                    min.cells = 0, min.features = 0)

    s2@meta.data$stim = s$stim[stim_perm]
    s2@meta.data$celltype = s$celltype[celltype_perm]
    s2 = process.dataset(s2, nfeatures = nrow(s2), npcs=k)

    # I need to recompute the snr because the library sizes have changed
    # so the normalized expression values have changed
    s2@misc$DE = s@misc$DE[is.element(s@misc$DE, rownames(s2))]
    s2@misc$snr = snr.analysis_util(s2, s2$stim, s2@misc$DE)

    saveRDS(s2, rds_file)
  }

  a = readRDS(rds_file)
  return (a)
}



nISG.dataset_DE = function(nISG, permutation_type, dataset, overwrite=F)
{
  rds_file = file.path(dir.dataset_DE(dataset),
                       paste(dataset, "_", nISG,
                             "_", permutation_type, ".rds", sep=""))

  if (!file.exists(rds_file) | overwrite) {

    s = permutation.dataset_DE(permutation_type, dataset)
    # copying code from Kang, so use ISG for a DE gene
    ISGs = s@misc$DE
    non_ISGs = setdiff(rownames(s), ISGs)

    if (nISG > 0)
      s2 = s[c(non_ISGs,ISGs[1:nISG]),]
    else
      s2 = s[non_ISGs,]

    s2 = process.dataset(s2, npcs=20)

    # I need to recompute the snr because the library sizes have changed
    # so the normalized expression values have changed
    s2@misc$DE = s@misc$DE[is.element(s@misc$DE, rownames(s2))]
    s2@misc$snr = snr.analysis_util(s2, s2$stim, s2@misc$DE)

    # get the svd of the baseline
    Y = scaled_counts.dataset(s2)
    uv = irlba::svdr(Y, k=40)
    s2@misc$uv = uv

    saveRDS(s2, rds_file)
  }

  s = readRDS(rds_file)
  return (s)

}


########################################################
# utility functions
gene_table.dataset_DE = function(dataset, overwrite=F)
{
  min_pval=0.05
  t_file = here::here(dir.dataset_DE(dataset), paste(dataset, "_gene_table.csv", sep=""))

  if (!file.exists(t_file) | overwrite) {
    s = get_unprocessed.dataset_DE(dataset)
    # if s is too big, subsample
    if (ncol(s) > 5000)
      s = s[,sample.int(ncol(s), 5000)]
    cellnames = colnames(s)
    c.stim = cellnames[s$stim=="stim"]
    c.ctrl = cellnames[s$stim=="ctrl"]

    d = Seurat::FindMarkers(s@assays$RNA,
                            cells.1=c.stim,
                            cells.2=c.ctrl,
                            min.pct=0.05, logfc.threshold=0) %>%
      dplyr::filter(p_val < min_pval)
    d = dplyr::mutate(d, gene=rownames(d), .before=1)
    write.csv(d, t_file, row.names = F)
  }

  d = read.csv(t_file, check.names = F)
  return (d)
}

