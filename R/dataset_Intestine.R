# Alexandria project dataset
# https://singlecell.broadinstitute.org/single_cell/study/SCP1547/murine-small-intestinal-organoid-xpo1-inhibition?cluster=UMAP&spatialGroups=--&annotation=cell_subset--group--cluster&subsample=all#study-summary
# downloaded raw count file mu_org_KPT_raw_expression_counts.csv.gz
# meta data file :mu_intest_metadata.csv


dir.Intestine = function()
{
  root_dir = here::here("data", "Intestine_dataset")

  return (root_dir)
}

unprocessed.Intestine <- function(overwrite=F)
{
  raw_counts = here::here(dir.Intestine(),
                          "mu_org_KPT_raw_expression_counts.csv")
  ann_f = here::here(dir.Intestine(), "mu_intest_metadata.csv")
  rds_f = here::here(dir.Intestine(), "unprocessed_Intestine.rds")

  if (!file.exists(rds_f) | overwrite) {

    ann = data.table::fread(ann_f)
    # the first row seems to give the column class
    ann = ann[2:nrow(ann),]
    D = data.table::fread(raw_counts)
    genes = D[,1]
    D = D[,-1] %>% as.matrix

    # check
    if (any(ann$NAME != colnames(D)))
      browser()

    celltypes = ann$cell_subset
    rownames(D) = genes
    a <- Seurat::CreateSeuratObject(counts = D,
                                    project = "Intestine",
                                    min.cells = 0, min.features = 0)
    a$celltype = celltypes
    a[["percent.mt"]] <- PercentageFeatureSet(a, pattern = "^mt-")
    a = subset(a, subset = nCount_RNA > 2000 & nCount_RNA < 15000 & percent.mt < 20)

    a = process.dataset(a)
    a[["stim"]] = rep("ctrl", ncol(a))
    a$stim[a$celltype == "enteroendocrine"] = "stim"

    saveRDS(a, rds_f)
  }

  a = readRDS(rds_f)
  return (a)
}


