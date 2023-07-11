#' Data from Ordavas-Montanas, Shalek
#' "“Allergic inflammatory memory in human respiratory epithelial progenitor cells“
#'
#' downloaded data from  https://shaleklab.com/resource/mapping-allergic-inflammation/
dir.HumanNasal = function()
{
  root_dir = here::here("data", "HumanNasal_dataset")

  return (root_dir)
}


#' Construct and load the Pezoldt Seurat objects
#'
#' @param overwrite if Seurat object has already been constructed, should
#' it be overwritten?
unprocessed.HumanNasal <- function(overwrite=F)
{
  raw_counts = here::here(dir.HumanNasal(), "20180822_PolypAll_cleaned_rawdata.txt")
  ann_f = here::here(dir.HumanNasal(), "meta.txt")
  rds_f = here::here(dir.HumanNasal(), "unprocessed_HumanNasal.rds")

  if (!file.exists(rds_f) | overwrite) {

    ann = data.table::fread(ann_f)
    D = data.table::fread(raw_counts)
    genes = D[,1]
    D = D[,-1] %>% as.matrix

    # check
    if (any(ann$V1 != colnames(D)))
      browser()

    celltypes = ann$subset
    rownames(D) = genes
    a <- Seurat::CreateSeuratObject(counts = D,
                                    project = "HumanNasal",
                                    min.cells = 0, min.features = 0)
    a$celltype = celltypes
    a[["percent.mt"]] <- PercentageFeatureSet(a, pattern = "^MT-")
    a = subset(a, subset = nCount_RNA > 2000 & nCount_RNA < 6000 & percent.mt < 5)

    a = process.dataset(a)
    a[["stim"]] = rep("ctrl", ncol(a))
    a$stim[a$celltype == "Ciliated"] = "stim"

    saveRDS(a, rds_f)
  }

  a = readRDS(rds_f)
  return (a)
}

