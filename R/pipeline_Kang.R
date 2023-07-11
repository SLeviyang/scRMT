dir.Kang_IFN = function(celltype=NULL, insilico=T)
{
  root_dir = here::here("analysis", "Kang_IFN_datasets")
  if (!dir.exists(root_dir))
    dir.create(root_dir)

  if (is.null(celltype))
    return (root_dir)

  if (insilico)
    celltype_dir = here::here("analysis", "Kang_IFN_datasets",
                            paste("insilico_", celltype, sep=""))
  else
    celltype_dir = here::here("analysis", "Kang_IFN_datasets",
                              celltype)
  if (!dir.exists(celltype_dir))
    dir.create(celltype_dir)

  return (celltype_dir)
}

dir_analysis.Kang_IFN = function()
{
  d = here::here("analysis", "Kang_IFN_tables")
  if (!dir.exists(d))
    dir.create(d)
  
  return (d)
}


filename.Kang_IFN = function(insilico, celltype, p, ISG_type=NULL,
                             permutation_type=NULL,
                             nISG=NULL)
{
  if (is.null(ISG_type))
    ISG_type = "cellspecific"

  f_name = ifelse(insilico,
                  paste("insilico_", celltype, "_", p, "_", ISG_type, ".rds",sep=""),
                  paste("Kang_", celltype, "_", p, "_", ISG_type, ".rds",sep=""))
  if (!is.null(nISG))
    f_name = paste("nISG_", nISG, "_", permutation_type, "_", f_name, sep="")
  else if (!is.null(permutation_type))
    f_name = paste("permutation_", permutation_type, "_", f_name, sep="")
  else
    f_name = paste("base_", f_name, sep="")

  rds_file = here::here(dir.Kang_IFN(celltype, insilico), f_name)

  return (rds_file)
}

get_spike_k.Kang_IFN = function(celltype)
{
  if (celltype == "Kang_all") {
    spike_k = 20
  } else if (celltype == "splatter_all") {
    spike_k = 6
  } else {
    if (celltype %in% get_celltypes.Kang()) {
      spike_k = 20
    } else if (celltype %in% get_celltypes.Kang_insilico()) {
      spike_k = 4
    }
    else
      stop("celltype not identified")
  }
  return (spike_k)
}


###############################

get_generated_table.Kang_IFN = function()
{
  md = rbind(data.frame(celltype=get_celltypes.Kang(), insilico=F,
                        stringsAsFactors = F),
             data.frame(celltype=get_celltypes.Kang_insilico(),
                        insilico=T, stringsAsFactors = F))

  d_all = plyr::adply(md, 1, function(d) {
    celltype = d$celltype
    insilico = d$insilico

    cell_dir = dir.Kang_IFN(celltype, insilico)

    if (!dir.exists(cell_dir))
      return (NULL)

    f = dir(cell_dir)
    f = f[grepl("nISG", f)]
    if (length(f)==0)
      return (NULL)

    # remove .rds
    f = sapply(f, function(s) substr(s, 1, nchar(s)-4))
    # I have problems with Kang_all and splatter_all
    f =  stringr::str_replace(f, pattern="_all", replacement="")

    fs = strsplit(f, split="_", fixed=T)

    nISG = sapply(fs, "[", 2) %>% as.numeric
    permutation_type = sapply(fs, "[", 3)
    insilico = sapply(fs, "[", 4) == "insilico"
    celltype = sapply(fs, "[", 5)
    p = sapply(fs, "[", 6) %>% as.numeric
    ISG_type = sapply(fs, "[", 7)

    d = data.frame(celltype=celltype,
                   insilico=insilico,
                   p=p,
                   permutation_type=permutation_type,
                   ISG_type=ISG_type,
                   nISG=nISG,stringsAsFactors = F) %>%
        dplyr::arrange(celltype, insilico, p,
                       permutation_type, ISG_type, nISG)
    rownames(d) = NULL
    return (d)
  })
  d_all$celltype = stringr::str_replace(d_all$celltype, "Kang", "Kang_all")
  d_all$celltype = stringr::str_replace(d_all$celltype, "splatter", "splatter_all")

  return (d_all)
}

make.Kang_IFN = function()
{
  make_both_celltype.Kang_IFN()
  make_all.Kang_IFN()
}

#' Make the "all" celltypes versions of insilico and Kang datasets
make_all.Kang_IFN = function()
{
  for (insilico in c(T,F)) {
    if (insilico)
      all_name = "splatter_all"
    else
      all_name = "Kang_all"
    s = base.Kang_IFN(insilico = insilico, all_name, 0.5, "joint")
    parallel::mclapply(c("none"), function(pt) {
      permutation.Kang_IFN(pt, insilico=insilico, all_name, 0.5, "joint")
      return (NULL)
    }, mc.cores=3)

    total_nISG = s@misc$ISG %>% length
    if (insilico)
      nISG_vec = seq(0, total_nISG, length.out=40) %>% round
    else
      nISG_vec = 0:total_nISG
    parallel::mclapply(nISG_vec, function(nISG) {
      cat("creating nISG", nISG, "\n")
    #  nISG.Kang_IFN(nISG, "orthogonal", insilico=insilico,
    #                all_name, 0.5, "joint")
      nISG.Kang_IFN(nISG, "none",  insilico=insilico,
                    all_name, 0.5, "joint")
    #  nISG.Kang_IFN(nISG, "single",  insilico=insilico,
    #                all_name, 0.5, "joint")
      return (NULL)
    }, mc.cores=3)
  }
}

make_both_celltype.Kang_IFN = function()
{
  make_celltype.Kang_IFN(T)
  make_celltype.Kang_IFN(F)

}

make_celltype.Kang_IFN = function(insilico)
{
  if (insilico)
      ct = get_celltypes.Kang_insilico() %>% setdiff(c("splatter_all", "ct5"))
  else
      ct = get_celltypes.Kang() %>% setdiff(c("Kang_all", "CD4 T cels"))

  pt = c("none")
  p = c(0.05, 0.25, 0.5)
  d = expand.grid(celltype=ct, p=p, permutation_type=pt)

  for (i in 1:nrow(d))
    make_nISG.Kang_IFN(d$permutation_type[i],
                       insilico=insilico,
                       celltype=d$celltype[i],
                       p=d$p[i],
                       ISG_type=NULL)

}

make_nISG.Kang_IFN = function(permutation_type, insilico, celltype,
                              p, ISG_type=NULL,
                              N=50, max_signal=6)
{
  bf = filename.Kang_IFN(insilico, celltype, p, ISG_type)
  if (!file.exists(bf))
    base.Kang_IFN(insilico, celltype, p, ISG_type)

  pf = filename.Kang_IFN(insilico, celltype, p, ISG_type, permutation_type)
  s = permutation.Kang_IFN(permutation_type, insilico, celltype, p, ISG_type)

  snr = snr.analysis_util(s, s$stim, s@misc$ISG)
  signal = cumsum(snr/(1+snr))

  max_nISG = which(signal <= max_signal) %>% max
  nISG = seq(0, max_nISG, length.out=N) %>% round %>% unique
  nISG = c(nISG, 0:min(20, max_nISG)) %>% unique %>% sort

  parallel::mclapply(nISG, function(i) {
    f = filename.Kang_IFN(insilico, celltype, p, ISG_type, permutation_type, i)

    if (!file.exists(f))
      s = nISG.Kang_IFN(i, permutation_type, insilico,
                        celltype, p, ISG_type)
    return (NULL)
  }, mc.cores=5)
}


########################################################
#' Processes the raw dataset for a (celltype, p) combination
#' Prior to processing, ISGs are restricted by type
#'
#' @param celltype
#' @param p one of 0.05, 0.25, 0.5
#' @param ISG_type  Only relevant if celltype is all.  In that case
#' can be "joint" or a celltype.  If a celltype then celltype specific ISGs are used
#'
base.Kang_IFN = function(insilico, celltype, p, ISG_type=NULL, overwrite=F)
{
  if (insilico)
    all_ct = get_celltypes.Kang_insilico() %>% setdiff("splatter_all")
  else
    all_ct = get_celltypes.Kang() %>% setdiff("Kang_all")

  if (celltype == "Kang_all" | celltype=="splatter_all") {
    if (is.null(ISG_type)) {
      message("If the celltype is all, then ISG_type must be specified")
      return (NULL)
    } else if (!is.element(ISG_type, c("joint", all_ct))) {
    message("Celltype is all, but ISG_type is not valid")
    return (NULL)
    }
  }

  if (celltype != "Kang_all" & celltype != "splatter_all" & !is.null(ISG_type)) {
    message("If the celltype is not all, then ISG_type should not be specified")
    return (NULL)
  }


  if (is.null(ISG_type))
    ISG_type = "cellspecific"

  rds_file = filename.Kang_IFN(insilico, celltype, p, ISG_type)
  if (!file.exists(rds_file) | overwrite) {

    if (insilico) {
      s = stim_freq.Kang_insilico(celltype, p)
      ISG_table = ISG_table.Kang_insilico(celltype)
    }
    else {
      s = stim_freq.Kang(celltype, p)
      ISG_table = ISG_table.Kang(celltype)
    }

    all_ISGs = intersect(rownames(s), ISG_table$gene)
    non_ISGs = setdiff(rownames(s), all_ISGs)

    if (ISG_type == "joint")
      active_ISGs = dplyr::filter(ISG_table, type=="joint") %>% dplyr::pull(gene)
    else if (ISG_type == "cellspecific") {
      ct = celltype
      active_ISGs = dplyr::filter(ISG_table, type=="cellspecific", celltype==ct) %>%
                    dplyr::pull(gene)
    }
    else
      stop("Problem!  We should never get here")

    s = s[c(non_ISGs, active_ISGs),]
    s = process.dataset(s, nfeatures=1000, scale.max=Inf, npcs=20)
    # randomize the ISG order for downstream analysis that is not
    # biased
    active_ISGs = intersect(rownames(s), active_ISGs)
    s@misc$ISG = active_ISGs[sample.int(length(active_ISGs))]
    # for joint ISG, order from least snr to greatest to resolve threshold

    if (ISG_type == "joint") {
      snr = snr.analysis_util(s, s$stim, s@misc$ISG)
      ind = order(snr, decreasing=F)
      s@misc$ISG = s@misc$ISG[ind]
    }

    # I don't need the spike to compute snr and Xs is computed
    # on the datasets rather than Kang_IFN
    if (insilico) {
      s = update_Seurat_splatter_info(s)
          #recompute_spike_from_Seurat()
    }

    s@misc$snr = snr.analysis_util(s, s$stim, s@misc$ISG)
    Y = scaled_counts.dataset(s)
    uv = irlba::svdr(Y[,s@misc$ISG], k=3)
    rownames(uv$v) = s@misc$ISG
    s@misc$ISG_uv = uv

    saveRDS(s, rds_file)
  }

  a = readRDS(rds_file)
  return (a)
}

#########################################################
# permutation of Kang celltype datasets

#' A Kang cell type dataset with permutations to control modules
#'
#' @param insilico
#' @param celltype
#' @param p frequency of stimulated cells
#' @param ISG_type
#' @param permutation_type
#'
#' @details After permutation of cells the celltype and stim labels
#' have a new meaning.   The stim label preserves cells that are stimulated
#' vs unstimulated.   The celltype label preserves reflects the gene
#' expression pattern of non-ISGs.   In "single", the celltype labels are
#' completely randomized and have no meaning.   In "orthogonal", the cell
#' type labels reflect shared non-ISG expression
permutation.Kang_IFN = function(permutation_type,
                                insilico, celltype, p, ISG_type=NULL,
                                k=20,
                                overwrite=F)
{
  if (!is.element(permutation_type, c("single", "orthogonal", "none"))) {
    message("invalid permutation type")
    return (NULL)
  }

  rds_file = filename.Kang_IFN(insilico, celltype, p, ISG_type, permutation_type)

  if (!file.exists(rds_file) | overwrite) {
    s = base.Kang_IFN(insilico, celltype, p, ISG_type)

    ISGs = s@misc$ISG

    DE = is.element(rownames(s), ISGs)

    ncells = ncol(s)

    if (permutation_type == "single") {
      per = permutation_single.Kang_IFN(DE, s$stim)
    }
    else if (permutation_type == "orthogonal") {
      per = permutation_ortho.Kang_IFN(ncells, DE)
    }
    else
      per = permutation_none.Kang_IFN(ncells, DE)

    # The permutations change the celltype,stim labels
    stim_perm = per[which(DE)[1],]
    celltype_perm = per[which(!DE)[1],]

    Y = s@assays$RNA@counts %>% Matrix::as.matrix() %>%
      apply_permutation.Kang_IFN(per)
    if (insilico) {
      X = s@misc$X %>% apply_permutation.Kang_IFN(per)
      s2 = Seurat_from_splatter_info(X, Y, s$celltype[celltype_perm],
                                     s@misc$DE, s@misc$gene_patterns,
                                     s@misc$dispersion,
                                     gene_names=rownames(s),
                                     cell_names=colnames(s),
                                     compute_Xs = F)


    }
    else
      s2 = Seurat::CreateSeuratObject(Y,
                                      project = "Kang_et_al",
                                      min.cells = 0, min.features = 0)

    s2@meta.data$stim = s$stim[stim_perm]
    s2@meta.data$celltype = s$celltype[celltype_perm]

    s2 = process.dataset(s2, nfeatures = nrow(s2), npcs=k)
    s2@misc$ISG = s@misc$ISG[is.element(s@misc$ISG, rownames(s2))]
    # I have to recompute the snr because library sizes have changes
    s2@misc$snr = snr.analysis_util(s2, s2$stim, s2@misc$ISG)

    s2@misc$ISG_uv = s@misc$ISG_uv
    s2@misc$ISG_uv$v = s2@misc$ISG_uv$v[s2@misc$ISG,]

    saveRDS(s2, rds_file)
  }

  a = readRDS(rds_file)
  return (a)
}

#' Constructs permutation matrix
#' DE genes are permuted within cell types
#' non DE genes are permuted without respect to cell type
#'
#' @param DE a logical vector giving the DE genes
#' @param celltypes a character vector giving the cell types
permutation_single.Kang_IFN = function(DE, celltypes)
{
  u_ct = unique(celltypes)
  ct_list = lapply(u_ct, function(cct) which(celltypes==cct)) %>%
            setNames(u_ct)
  ncells = length(celltypes)

  p = sapply(1:length(DE), function(i) {
    g = 1:ncells

    if (DE[i])
      for (cct in u_ct) {
        ind = ct_list[[cct]]
        g[ind] = g[sample(ind)]
      }
    else
      g = g[sample.int(ncells)]

    return (g)
  }) %>% t

  return (p)
}


#' Permutes genes jointly to remove correlations between the DE and non-DE
#' genes.  non DE genes are permuted jointly without respect to cell type
#'
#' @param DE a logical vector giving the DE genes
permutation_ortho.Kang_IFN = function(ncells, DE)
{
  sample_p = sample.int(ncells)
  g = 1:ncells

  p = sapply(1:length(DE), function(i) {
    if (!DE[i])
     return (sample_p)
    else
      return (g)
  }) %>% t

  return (p)
}

#' No permutation
#'
#' @param DE a logical vector giving the DE genes
permutation_none.Kang_IFN = function(ncells, DE)
{
  ngenes = length(DE)
  g = 1:ncells
  p = replicate(ngenes, g) %>% t

  return (p)
}


apply_permutation.Kang_IFN = function(X, p)
{
  if (any(dim(X) != dim(p)))
    stop("problem")

  for (i in 1:nrow(X))
    X[i,] = X[i,p[i,]]

  return (X)
}


#####################################################################
# subset frequency of stimulated cells

#' Subset the ISGs of a particular Kang dataset
#'
#' @param celltype
#' @param p frequency of stimulated cells
#' @param permutation_type "single", "orthogonal" or "type"
#' @param nISG an integer
#' @param in_silico
#' @param k number of pcs
nISG.Kang_IFN = function(nISG,
                         permutation_type,
                         insilico,
                         celltype,
                         p,
                         ISG_type=NULL,
                         k=20, overwrite=F)
{
  if (!is.element(permutation_type, c("single", "orthogonal", "none"))) {
    message("invalid permutation type")
    return (NULL)
  }

  rds_file = filename.Kang_IFN(insilico, celltype, p, ISG_type, permutation_type, nISG)

  if (!file.exists(rds_file) | overwrite) {
    s = permutation.Kang_IFN(permutation_type, insilico, celltype,
                             p, ISG_type)
    ISGs = s@misc$ISG
    non_ISGs = setdiff(rownames(s), ISGs)

    if (nISG > 0)
      s2 = s[c(non_ISGs,ISGs[1:nISG]),]
    else
      s2 = s[non_ISGs,]

    s2 = process.dataset(s2, npcs=k)
    if (insilico & nISG==0) {
      s2 = update_Seurat_splatter_info(s2) %>%
           recompute_spike_from_Seurat()
    }

    s2@misc$ISG = intersect(rownames(s2), s@misc$ISG)
    # I have to recompute the snr because library sizes have changes
    s2@misc$snr = snr.analysis_util(s2, s2$stim, s2@misc$ISG)

    s2@misc$ISG_uv = s@misc$ISG_uv
    s2@misc$ISG_uv$v = s2@misc$ISG_uv$v[s2@misc$ISG,]

    Y = scaled_counts.dataset(s2)
    uv = irlba::svdr(Y, k=k)
    s2@misc$uv = uv
    saveRDS(s2, rds_file)
  }

  a = readRDS(rds_file)
  return (a)
}

