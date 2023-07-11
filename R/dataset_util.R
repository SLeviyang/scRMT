#' Processing of Seurat dataset
#' Applies standard Seurat processing to a dataset
#'
#' @param a Seurat object with non-empty assay$data
#' @param nfeatures number of genes to include
#' @param k number of nearest neighbors to use in umap
#' @param npcs number of principal components.
#' @param resolution resolution for Seurat cell clustering algorithm
#' @param scale.max parameter to Seurat ScaleData function
#'
#' @returns A Seurat object with genes limited to variable genes.
#' Seurat object will include scaled data, pca and umap reductions,
#' and clustering.
process.dataset = function(a, nfeatures=1000, k=10, npcs=20,
                          resolution=0.5,
                          scale.max=10)
{
  a = Seurat::FindVariableFeatures(a, selection.method = "vst",
                                   nfeatures = min(nrow(a), nfeatures))
  a = a[Seurat::VariableFeatures(object = a),]

  a = Seurat::NormalizeData(a,
                            normalization.method="LogNormalize",
                            scale.factor=10000)
  a = Seurat::ScaleData(a, features = rownames(a), scale.max=scale.max)

  a = Seurat::RunPCA(a,
                     features = rownames(a),
                     npcs=npcs)

  a <- Seurat::FindNeighbors(a, dims = 1:npcs)
  a <- Seurat::FindClusters(a, resolution = resolution,
                            graph.name="RNA_nn")
  a = Seurat::RunUMAP(a, dims=1:npcs)
  # build a knn
  #X = a@reductions$pca@cell.embeddings
  #K =  BiocNeighbors::findKNN(X, k=k, BNPARAM=BiocNeighbors::KmknnParam())
  #a@misc$K = K

  return (a)
}

#' Returns an expression matrix with genes centered and scaled to norm 1.
#'
#' @param a Seurat object
#' @param use.scale.data should the Seurat RNA@data matrix be scaled
#' or should RNA@scale.data
#'
#' @return a scaled matrix with genes as columns
#'
#' @details Seurat's scale.data matrix is not actually scaled because
#' expression counts above 10 are rounded down to 10 (see scale.max
#' parameter in Seurat::ScaleData).  Also, genes are scaled to unit variance,
#' not unit norm in Seurat, although there too the rounding means the
#' columns are not actually scaled.
scaled_counts.dataset = function(a,
                                 use.scale.data=T)
{
  if (use.scale.data)
    counts = a@assays$RNA@scale.data %>% t
  else
    counts = as.matrix(a@assays$RNA@data) %>% t
  counts = base::scale(counts, center=T)
  counts = apply(counts, 2, function(g) g/sqrt(sum(g^2)))

  rownames(counts) = colnames(a)
  colnames(counts) = rownames(a)

  return (counts)
}

plot_feature.dataset = function(a, v,
                                reduction="umap",
                                comparison=NULL)
{
  if (!is.null(comparison))
    par(mfrow=c(1,2))
  a@meta.data$feature = v
  g1 = Seurat::FeaturePlot(a, "feature", reduction=reduction)

  if (!is.null(comparison))
    g2 = Seurat::DimPlot(a, group.by=comparison)

  if (!is.null(comparison))
    gridExtra::grid.arrange(g1, g2, ncol=2)
  else
    gridExtra::grid.arrange(g1, ncol=1)

  return (NULL)
}


####################################################################
# permutation functions


#' Permutes the DE genes within cell type.  The non DE genes
#' are permuted independently so the notion of a cell type is lost.
permutation_single.dataset_util = function(DE, celltypes)
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



#' Constructs permutation matrix
#' DE genes are permuted within cell types
#' non DE genes are not changed (this differs from Kang)
#'
#' @param DE a logical vector giving the DE genes
#' @param celltypes a character vector giving the cell types
permutation_orthogonal.dataset_util = function(ncells, DE)
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
permutation_none.dataset_util = function(ncells, DE)
{
  ngenes = length(DE)
  g = 1:ncells
  p = replicate(ngenes, g) %>% t

  return (p)
}



#' Don't change baseline,  but resample perturbation module
permutation_module.dataset_util = function(DE, celltypes)
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

    return (g)
  }) %>% t

  return (p)
}
