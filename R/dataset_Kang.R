# Kang HM, Subramaniam M, Targ S, Nguyen M et al.
#Multiplexed droplet single-cell RNA-sequencing using
#     natural genetic variation. Nat Biotechnol 2018
# accession: GSE96583
#
# Run 2 control:  pooled PBMC from 8 lupus patients,
#                 untreated for 6 hours (GSM2560248)
# Run 2 IFNB:  pooled PBMC from 8 lupus patients,
#               treated with IFN B for 6 hours (GSM2560249)

dir.Kang = function()
{
  root_dir = here::here("data", "Kang_dataset")
  if (!dir.exists(root_dir))
    dir.create(root_dir)

  return (root_dir)
}

load_annotations.Kang = function()
{
  root_dir = here::here(dir.Kang(), "Kang_raw_files")
  annotation_file = file.path(root_dir, "annotation.tsv")

  return (read.csv(annotation_file, sep="\t"))
}

make.Kang = function()
{
  make_base_datasets.Kang()
  make_all.Kang()
  make_celltype.Kang()
}

#' Processes GEO files for the Kang dataset
#'
#' @details starts with GEO files in Kang_GEO_dataset folder
make_base_datasets.Kang = function()
{
  GEO_conversion.Kang()
  unprocessed.Kang()
}

make_all.Kang = function()
{
  celltype.Kang("Kang_all")
  processed_celltype.Kang("Kang_all")
  stim_freq.Kang("Kang_all", 0.5)
  ISG_table.Kang("Kang_all")
}

make_celltype.Kang = function()
{
  ct = get_celltypes.Kang() %>% setdiff(c("Kang_all", "CD4 T cells"))
  parallel::mclapply(ct, function(celltype) {
    cat("making cell type file", celltype, "\n")
    rds_file = here::here(dir.Kang(),
                          paste("Kang_", celltype, ".rds",sep=""))


    if (!file.exists(rds_file))
      celltype.Kang(celltype)
    return (NULL)
  },mc.cores = 5)

  parallel::mclapply(ct, function(celltype) {
    cat("making processed cell type file", celltype, "\n")
    rds_file = here::here(dir.Kang(),
                          paste("processed_Kang_", celltype, ".rds",sep=""))


    if (!file.exists(rds_file))
      processed_celltype.Kang(celltype)
    return (NULL)
  }, mc.cores=5)

  p = c(0.05, 0.25, 0.5)
  d = expand.grid(celltype=ct, p=p, stringsAsFactors = F)

  parallel::mclapply(1:nrow(d), function(i) {
    print(d[i,])

    celltype = d$celltype[i]
    p = d$p[i]
    rds_file = here::here(dir.Kang(),
                          paste("Kang_", celltype, "_", p, ".rds",sep=""))

    if (!file.exists(rds_file))
      stim_freq.Kang(celltype, p)
    return (NULL)
  },mc.cores = 5)

  for (cct in ct)
    ISG_table.Kang(cct)

  return (NULL)
}

get_celltypes.Kang = function(full=F, overwrite=F)
{
  tb_file = here::here(dir.Kang(), "cell_types.csv")

  if (!file.exists(tb_file) | overwrite) {
    fs = unprocessed.Kang()
    d = fs$cell %>% table %>% as.data.frame %>%
      setNames(c("celltype", "n"))
    d = rbind(d, data.frame(celltype="Kang_all", n=sum(d$n)))

    write.csv(d, tb_file, row.names = F)
  }

  ct = read.csv(tb_file)
  if (full)
    return (ct)
  else
    return (ct$celltype)
}


ISG_table.Kang = function(celltype, overwrite=F)
{
  min_pval=0.05
  t_file = here::here(dir.Kang(),
                      paste("Kang_ISG_table_", celltype, ".csv", sep=""))

  if (!file.exists(t_file) | overwrite) {

    if (celltype != "Kang_all") {
      # consider ISG from top 2000 to avoid low expression counts that
      # may mess up DE computations
      s = celltype.Kang(celltype) %>% process.dataset(nfeatures=2000)
      cellnames = colnames(s)
      cells.stim = cellnames[s$celltype == celltype & s$stim=="stim"]
      cells.ctrl = cellnames[s$celltype == celltype & s$stim=="ctrl"]
      d = Seurat::FindMarkers(s@assays$RNA,
                                cells.1=cells.stim,
                                cells.2=cells.ctrl,
                                min.pct=0.05, logfc.threshold=0) %>%
          dplyr::filter(p_val < min_pval)
      d = dplyr::mutate(d, type="cellspecific",
                        celltype=celltype, gene=rownames(d), .before=1)
    } else {
      # I assume that the Kang datasets restricted to a single cell type
      # have been computed, and I use those to create a ISG table for
      # all cell types
      ct = get_celltypes.Kang() %>% setdiff("Kang_all")
      n_ct = length(ct)
      d_in = plyr::adply(ct, 1, ISG_table.Kang, .id = NULL)
      tg = table(d_in$gene)

      cellspecific_genes = which(tg == 1) %>% names
      joint_genes = which(tg == n_ct) %>% names
      other_genes = setdiff(d_in$gene, c(joint_genes, cellspecific_genes))
      d_cs = d_in[is.element(d_in$gene, cellspecific_genes),]
      d_j = data.frame(type="joint", celltype=NA,
                       gene=joint_genes, p_val=NA, avg_log2FC=NA,
                       pct.1=NA, pct.2=NA, p_val_adj=NA, stringsAsFactors = F)
      d_other = data.frame(type="other", celltype=NA,
                           gene=other_genes, p_val=NA, avg_log2FC=NA,
                           pct.1=NA, pct.2=NA, p_val_adj=NA, stringsAsFactors = F)
      d = rbind(d_j, d_cs, d_other)
    }

    write.csv(d, t_file, row.names = F)
  }


  d = read.csv(t_file, header=T, check.names = F)
  return (d)
}

###########################################################################
#' Filters simulated cells to make the stimulated
#' cells a frequency p of the total number of cells
#'
#' @param celltype
#' @param p num_stim/total
stim_freq.Kang = function(celltype, p, overwrite=F)
{
  if (p > 0.5) {
    message("p must not exceed 0.5")
    return (NULL)
  }

  rds_file = here::here(dir.Kang(),
                        paste("Kang_", celltype, "_", p, ".rds",sep=""))

  if (!file.exists(rds_file) | overwrite) {
    s = celltype.Kang(celltype)
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

    saveRDS(s, rds_file)
  }

  a = readRDS(rds_file)
  return (a)
}

#' Performs Seurat workflow on the celltype.Kang output.  Also
#' computes an svd for downstream analysis.
processed_celltype.Kang = function(celltype, overwrite=F)
{
  rds_file = here::here(dir.Kang(),
                        paste("processed_Kang_", celltype, ".rds",sep=""))

  if (!file.exists(rds_file) | overwrite) {
    s = celltype.Kang(celltype) %>%
        process.dataset(nfeatures=1000, scale.max=Inf, npcs=20)
    Y = scaled_counts.dataset(s)
    s@misc$uv = irlba::svdr(Y, k=20)

    saveRDS(s, rds_file)
  }

  a = readRDS(rds_file)
  return (a)
}


#' Filters the full Kang dataset for a cell type.
#' The cells are filtered in a minimial way
#' to make the number of stimulated and control cells equal.
celltype.Kang = function(celltype, overwrite=F)
{
    rds_file = here::here(dir.Kang(),
                          paste("Kang_", celltype, ".rds",sep=""))

    if (!file.exists(rds_file) | overwrite) {
      s = unprocessed.Kang()
      if (celltype != "Kang_all")
        s = s[,s$cell == celltype]
      s$celltype = s$cell
      s@meta.data$cell = NULL

      # make all the cell types have equal stim and ctrl for
      # the benifit of downstream manipulations
      all_ct = s$celltype %>% unique

      ind_list = lapply(all_ct, function(ct) {
        ind_stim = which(s$stim == "stim" & s$celltype == ct)
        ind_ctrl = which(s$stim == "ctrl" & s$celltype == ct)
        if (length(ind_stim) > length(ind_ctrl))
          ind_stim = sample(ind_stim, length(ind_ctrl))
        else
          ind_ctrl = sample(ind_ctrl, length(ind_stim))

        return (list(stim=ind_stim, ctrl=ind_ctrl))
      })
      ind_stim = lapply(ind_list, "[[", "stim") %>% unlist
      ind_ctrl = lapply(ind_list, "[[", "ctrl") %>% unlist

      s = s[,c(ind_stim, ind_ctrl)]
      #s =  process.dataset(s, nfeatures=1000, npcs=20, scale.max=Inf)

      s = RenameCells(s, new.names=paste(s$celltype, 1:ncol(s), sep="-"))
      saveRDS(s, rds_file)
    }

    a = readRDS(rds_file)
    return (a)
}


#' Creates and saves a Seurat object for the full Kang dataset
#' Includes stimulated and unstimulated PBMC
unprocessed.Kang <- function(overwrite=F)
{
  root_dir = here::here(dir.Kang(), "Kang_raw_files")
  stim_dir = file.path(root_dir, "stimulated")
  control_dir = file.path(root_dir, "control")
  rds_f = file.path(dir.Kang(), "unprocessed_Kang.rds")

  if (file.exists(rds_f) & !overwrite)
    return (readRDS(rds_f))

  annotation = read.csv(file.path(root_dir, "annotation.tsv"),
                        sep="\t", header=T,
                        stringsAsFactors=F)
  annotation = dplyr::mutate(annotation,
                             barcode=strsplit(rownames(annotation), split="-") %>%
                               sapply("[", 1))
  stims = annotation$stim

  load_agent_dataset <- function(data_dir,
                                 stim_type) {
    bc = read.csv(file.path(data_dir,"barcodes.tsv"),
                  sep="\t", header=F, stringsAsFactors = F) %>%
      dplyr::pull("V1") %>%
      strsplit(rownames(annotation), split="-") %>%
      sapply("[", 1)

    d <- Seurat::Read10X(data.dir = data_dir)
    a <- Seurat::CreateSeuratObject(counts = d,
                                    project = "Kang_et_al",
                                    min.cells = 0, min.features = 0) %>%
      Seurat::AddMetaData(bc, col.name="barcode")

    md = dplyr::filter(annotation, stim == stim_type) %>%
      dplyr::inner_join(a@meta.data, by="barcode")
    if (nrow(md) != ncol(a))
      stop("problem")

    a@meta.data = md
    return (a)
  }

  a_stim = load_agent_dataset(stim_dir, "stim")
  a_ctrl = load_agent_dataset(control_dir, "ctrl")
  a = merge(a_ctrl, y=a_stim,
            project = "Kang")

  # remove multiples in droplets
  a = a[,a@meta.data$multiplets == "singlet"]
  #VISION: "cluster 3 (annotation from original study) were excluded
  # as these appeared to be proliferating T cells whose large difference
  # from the rest of the cells served to mask more fine-grained heterogeneity
  a = a[,a@meta.data$cluster != 3]

  a = a[,!is.na(a$cell)]
  a = a[,is.element(a$stim,c("stim", "ctrl"))]

  saveRDS(a, rds_f)
  return (a)
}


#' Organizes the raw GEO datafiles.  Included here for reproducibility
GEO_conversion.Kang <- function()
{
  raw_dir=here::here(dir.Kang(), "Kang_GEO_files")
  root_dir = here::here(dir.Kang(), "Kang_raw_files")

  if (!dir.exists(root_dir))
    dir.create(root_dir)

  # copy data from GEO filenames to informative file names
  file.copy(file.path(raw_dir,"GSE96583_batch2.total.tsne.df.tsv"),
            file.path(root_dir, "annotation.tsv"))
  annotation_file = file.path(root_dir, "annotation.tsv")

  stim_dir = here::here(root_dir, "stimulated")
  if (!dir.exists(stim_dir)) {
    dir.create(stim_dir)
    file.copy(file.path(raw_dir, "GSM2560249_2.2.mtx"),
              file.path(stim_dir, "matrix.mtx"))
    file.copy(file.path(raw_dir, "GSM2560249_barcodes.tsv"),
              file.path(stim_dir, "barcodes.tsv"))
    file.copy(file.path(raw_dir, "GSE96583_batch2.genes.tsv"),
              file.path(stim_dir, "genes.tsv"))
  }


  control_dir = here::here(root_dir, "control")
  if (!dir.exists(control_dir)) {
    dir.create(control_dir)
    file.copy(file.path(raw_dir, "GSM2560248_2.1.mtx"),
              file.path(control_dir, "matrix.mtx"))
    file.copy(file.path(raw_dir, "GSM2560248_barcodes.tsv"),
              file.path(control_dir, "barcodes.tsv"))
    file.copy(file.path(raw_dir, "GSE96583_batch2.genes.tsv"),
              file.path(control_dir, "genes.tsv"))
  }

  return (list(stimulated_dir=stim_dir,
               control_dir=control_dir,
               annotation_file=annotation_file))
}
