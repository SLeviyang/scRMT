# This is the 10x dataset used in the Seurat introductory vignette.
#Available at,
#`https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz`
#and with annotation information,
#`https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis`
#
# paper: Zheng et al., "Massively parallel digital transcriptional profiling of single cells".

dir.Zheng = function()
{
  root_dir = here::here("data", "Zheng_68K_PBMC_dataset")

  return (root_dir)
}



#' Construct and load the Zheng Seurat objects
#'
#' @param overwrite if Seurat object has already been constructed, should
#' it be overwritten?
unprocessed.Zheng <- function(overwrite=F, N=Inf)
{
  root_dir = dir.Zheng()
  annotation_file = file.path(root_dir, "68k_pbmc_barcodes_annotation.tsv")
  rds_f = file.path(root_dir, "unprocessed_Zheng.rds")

  if (!file.exists(rds_f) | overwrite) {
    set.seed(123)
    ann = read.table(annotation_file, sep="\t", header=T,
                     stringsAsFactors = F) %>%
          dplyr::select(barcodes, celltype)
    d <- Seurat::Read10X(data.dir = root_dir)
    a <- Seurat::CreateSeuratObject(counts = d,
                                    project = "Zheng_PBMC",
                                    min.cells = 0, min.features = 0)
    a@meta.data = cbind(a@meta.data, ann)

    # check
    if (!all(a@meta.data$barcodes == colnames(a)))
      stop("Zheng barcodes do not match annotation!")

    # subset for fast computation
    if (N < nrow(a))
      a = a[,sample.int(nrow(a), N)]

    # QC according to Seurat PBMC vignette
    a[["percent.mt"]] = Seurat::PercentageFeatureSet(a, pattern = "^MT-")
    a = subset(a, subset = nFeature_RNA > 200 &
                 nFeature_RNA < 2500 &
                 percent.mt < 5)

    # set up B cells as the module

    a = process.dataset(a)
    # the B cells defined in the annotation do not form a coherent cluster,
    # so here I restrict to cells which do.
    B_clust_ind = which(a$celltype=="CD19+ B" &
                          a@reductions$umap@cell.embeddings[,1] < -10 &
                          a@reductions$umap@cell.embeddings[,2] < 5)
    a[["stim"]] = rep("ctrl", ncol(a))
    a$stim[B_clust_ind] = "stim"

    saveRDS(a, rds_f)
  }

  a = readRDS(rds_f)
  return (a)
}
