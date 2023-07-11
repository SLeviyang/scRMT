
#' From a Seurat object using splatter simulation
#'
#' @param celltype_sizes number of cells in each cell type. If the vector
#' is named, then those names are used for the cell types, otherwise generic "ct#"
#' is used
#' @param patterns a list of logical vectors or a character vector.
#' If a list, each logical vector
#' has the same length as celltype_sizes and specifies whether a particular
#' cell type is differentially expressed (T) or at baseline (F).  If a character
#' vector, then entries must be binary, e.g. "0010",
#' with a 1 signifying DE and 0 baseline and string length equal to the number of
#' cell types.
#' @param ngenes a vector of integers of the same length as the patterns list.
#' Each entry gives the number of genes with the given DE pattern.
#' @param dispersion dispersion value for the NB distribution
#' @param single_DE_factor should all differential expression have the same
#' factor, which determines the mean up/down regulation?
#'
#' Returns a Seurat object.  The counts is in s@assays$RNA@counts.  Further
#' information can be found in s@misc through the following fields.
#' @field X mean expression for each gene by cell
#' @field DE differential expression factors (matrix of dimension genes by cell types)
#' @field gene_patterns a character vector showing the DE pattern used to create each gene
splatter = function(celltype_sizes,
                    patterns,
                    ngenes,
                    dispersion=10,
                    single_DE_factor=NULL,
                    de.scale=NULL,
                    de.mean=NULL,
                    debug=F)
{
  n_celltypes = length(celltype_sizes)

  if (length(ngenes) != length(patterns)) {
    message("patterns and ngenes must have the same length")
    return (NULL)
  }
  if (is.character(patterns)) {
    if (!all(nchar(patterns)==n_celltypes)) {
      message("invalid pattern strings")
      message("pattern length does not match number of cell types")
      if (debug) browser()
      return (NULL)
    }
    patterns = lapply(strsplit(patterns, split=""), function(p) p=="1")
  } else if (is.list(patterns)) {
    if (!all(sapply(patterns, length) == n_celltypes)) {
      message("invalid pattern logical vectors")
      message("pattern length does not match number of cell types")
      return (NULL)
    }
  } else {
    message("pattern must be a list or character vector")
    return (NULL)
  }

  ncells = sum(celltype_sizes)

  if (is.null(names(celltype_sizes)))
    names(celltype_sizes) = paste("ct", 1:n_celltypes, sep="")

  celltype_names = names(celltype_sizes)
  celltype = rep(celltype_names, celltype_sizes)


  s_list = mapply(function(p, ng) {
    names(p) = celltype_names
    splatter_module(celltype, ngenes=ng, de.pattern=p,
                    dispersion=dispersion,
                    single_DE_factor = single_DE_factor,
                    de.scale=de.scale,
                    de.mean=de.mean)
  }, patterns, ngenes, SIMPLIFY=F)

  X = do.call(rbind, lapply(s_list, "[[", "X"))
  Y = do.call(rbind, lapply(s_list, "[[", "Y"))
  DE = do.call(rbind, lapply(s_list, "[[", "DE"))

  s = Seurat_from_splatter_info(X, Y, celltype, DE,
                                gene_patterns_from_DE.splatter(DE),
                                dispersion)

  return (s)
}


#' Constructs splatter expressions for a gene module
#'
#' A gene module is an expression matrix in which all genes
#' have the same differential expression pattern across cell types
#'
#' @param celltype cell types of cells as a character.vector with
#' length equal to the total number of cells
#' @param ngenes number of genes to simulate
#' @param de.pattern a logical vector with a value for each
#' unique cell type, must be named using the cell types.  T means
#' de, F means baseline
#' @param dispersion negative binomial dispersion value
#' @param single_DE_factor should all genes have the identical
#' de factors, meaning identical means across the cells?  Not realistic,
#' but useful for theoretical studies.  If NULL then genes have
#' independent factors, otherwise should be a positive number which
#' is used as the factor
#'
#' returns a list
#' @field X mean matrix (genes by cells)
#' @field Y counts matrix (genes by cells)
#' @field DE differential expression factors (genes by cell types)
splatter_module = function(celltype,
                           ngenes,
                           de.pattern,
                           dispersion,
                           single_DE_factor=NULL,
                           de.mean=NULL,
                           de.scale=NULL)
{
  celltype_u = celltype %>% unique %>% sort
  ncell = length(celltype)
  if (length(celltype_u) != length(de.pattern)) {
    message("de.pattern does not have an entry for every cell type")
    return (NULL)
  }
  if (!setequal(celltype_u, names(de.pattern))) {
    message("names of de.pattern do not match the cell types")
    return (NULL)
  }
  if (!is.logical(de.pattern)) {
    message("de.pattern must be a logical vector")
    return (NULL)
  }


  de.pattern = de.pattern[celltype_u]
  # which cells should have DE
  de_cells = is.element(celltype, celltype_u[de.pattern]) %>% which
  non_de_cells = setdiff(1:ncell, de_cells)
  n_de_cells = length(de_cells)

  sce <- scater::mockSCE()
  group.prob = c(n_de_cells/ncell, 1 - n_de_cells/ncell)
  # splatter crashes if any group.prob is 0
  if (any(group.prob == 0)) {
    group.prob = group.prob + 1E-6
    group.prob = group.prob/sum(group.prob)
  }

  # make a few more cells than needed so we can get exact amount, splatter
  # assigns cells through bernoulli
  ncell_sim = round(1.25*ncell)
  params <- splatter::splatEstimate(sce)
  repeat {
    params <- splatter::setParams(params, update=list(nGenes=ngenes,
                                                    batchCells = ncell_sim,
                                                    out.prob=0,
                                                    dropout.type="none",
                                                    lib.loc=10000,
                                                    lib.scale=0,
                                                    lib.norm=T,
                                                    de.prob=c(1,0),
                                                    group.prob=group.prob
                                                  #  de.facLoc = .1*log(2)
                                                    #  de.downProb=0,
                                                    #  de.facScale = .02
    ))
    if (!is.null(de.scale))
      params <- splatter::setParams(params,list(de.facScale=de.scale))
    if (!is.null(de.mean))
      params <- splatter::setParams(params,list(de.facLoc=de.mean))

    sim <- splatter::splatSimulateGroups(params)

    # have we generated enough of each cell type?
    celltype_sim = sim@colData$Group %>% as.character
    de_cells_sim = which(celltype_sim == "Group1")
    non_de_cells_sim = which(celltype_sim == "Group2")

    if (length(de_cells_sim) >= length(de_cells) &
        length(non_de_cells_sim) >= length(non_de_cells))
      break

    ncell_sim = round(1.2*ncell_sim)
  }
  # BaseCellMeans have cell columns normalized to sum to lib_sizes
  X_sim = sim@assays@data$BaseCellMeans


  X = matrix(0, nrow=ngenes, ncol=ncell)
  X[,de_cells] = X_sim[,sample(de_cells_sim, length(de_cells))]
  X[,non_de_cells] = X_sim[,sample(non_de_cells_sim, length(non_de_cells))]

  # Y adds noise to X
  Y = matrix(rnbinom(nrow(X)*ncol(X),
                     size=dispersion,
                     mu=as.numeric(X)), ncol=ncol(X), nrow=nrow(X))


  de_factors = rowData(sim)[,"DEFacGroup1"]
  DE = matrix(1, nrow=ngenes, ncol=length(celltype_u))
  colnames(DE) = paste("DEFac", celltype_u, sep="")
  DE[,de.pattern] = de_factors

  gene_names = paste("gene", 1:nrow(X), sep="")
  cell_names = paste("cell", 1:ncell, sep="")
  rownames(X) = rownames(Y) = gene_names
  colnames(X) = colnames(Y) = cell_names
  rownames(DE) = gene_names

  return (list(X=X, Y=Y, DE=DE))
}


Seurat_from_splatter_info = function(X, Y, celltype, DE, gene_patterns,
                                     dispersion,
                                     gene_names=NULL,
                                     cell_names=NULL,
                                     debug=F,
                                     compute_Xs=T)
{
  if (is.null(gene_names))
    gene_names = paste("gene", 1:nrow(X), sep="")
  if (is.null(cell_names))
    cell_names = paste("cell", 1:ncol(X), sep="")
  rownames(X) = rownames(Y) = rownames(DE) = gene_names
  colnames(X) = colnames(Y) = cell_names
  names(gene_patterns) = gene_names
  names(celltype) = cell_names

  s <- CreateSeuratObject(counts = Y,
                          project = "splatter",
                          min.cells = 0, min.features = 0)
  s@meta.data$celltype = celltype
  s@misc$X = X
  s@misc$DE = DE
  s@misc$gene_patterns = gene_patterns
  s@misc$celltype = celltype
  s@misc$dispersion = dispersion

  # spike computations
  if (compute_Xs) {
    Xs = Seurat_spike(s@misc$X, dispersion=dispersion, debug=debug)
    rownames(Xs) = gene_names
    colnames(Xs) = cell_names
    em = rowSums(Xs^2) %>% sqrt
    s@misc$Xs = Xs
    s@misc$effective_mu = em
  } else
    s@misc$Xs = NULL

  return (s)
}

recompute_spike_from_Seurat = function(s)
{
  # the library sizes have changed so I need to recompute the spike
  Xs = Seurat_spike(s@misc$X, dispersion=s@misc$dispersion)

  rownames(Xs) = rownames(s)
  colnames(Xs) = colnames(s)

  s@misc$Xs = Xs
  s@misc$effective_mu = rowSums(Xs^2) %>% sqrt

  return (s)
}


#' Update the misc info when Seurat rows or columns have changed
update_Seurat_splatter_info = function(s)
{

  s@misc$X = s@misc$X[rownames(s), colnames(s)]
  s@misc$Xs = s@misc$Xs[rownames(s), colnames(s)]
  s@misc$DE = s@misc$DE[rownames(s),]
  s@misc$gene_patterns = s@misc$gene_patterns[rownames(s)]
  s@misc$cell_types = s@misc$cell_types[colnames(s)]

  s@misc$effective_mu = s@misc$effective_mu[rownames(s)]

  return (s)
}

###############################################
#' Merges genes over two splatter object assuming equal number of cells
#'
#' @param s1,s2 Seurat objects
splatter_rbind = function(s1, s2)
{

}
#' Split a cell type into two
splatter_split = function(s, target_group, q,
                          ngenes, dispersion=10)
{
  X = s@misc$X
  Y = s@assays$RNA@counts
  celltype = s@meta.data$celltype
  DE = s@misc$DE

  ngroups = ncol(DE)
  new_group = paste("Group", ngroups+1, sep="")
  # update DE so that the new group is identical to the target group
  tDE = DE[,paste("DEFac", target_group,sep="")]
  cnames = c(colnames(DE), paste("DEFac", new_group,sep=""))
  DE = cbind(DE,tDE)
  colnames(DE) = cnames

  # find fraction of target cell type
  p = mean(celltype == target_group)

  # in the split splatter object, Group1 is the new group and Group2
  #  are all other cell types
  s_split = splatter_DE(p=c(q*p, 1-q*p), de.prob=.5, dispersion=dispersion,
                      ncells=ncol(s), ngenes=ngenes)

  X_split = s_split@misc$X
  Y_split = s_split@assays$RNA@counts
  split_celltype = s_split@meta.data$celltype
  split_DE = s_split@misc$DE

  # put all of the split Group1 cell type into the target cell types rows, and
  # the Group2 cell types occupy the rest of the rows
  ind1 = which(split_celltype == "Group1")
  ind2 = setdiff(1:ncol(s), ind1)

  target_cell_ind = which(celltype == target_group)
  p_ind1 = target_cell_ind[1:length(ind1)]
  p_ind2 = setdiff(1:ncol(s), p_ind1)

  map_df = data.frame(split_index=c(ind1, ind2),
                      index=c(p_ind1, p_ind2),
                      split_celltype=split_celltype) %>%
           dplyr::arrange(index)
  map_ind = map_df %>% dplyr::pull(split_index)
  new_celltype_ind = which(map_df$split_celltype=="Group1")


  X_new = X_split[,map_ind]
  Y_new = Y_split[,map_ind]

  DE_t = split_DE[,"DEFacGroup2"]
  DE_add = matrix(DE_t, nrow=ngenes, ncol=ncol(DE))
  colnames(DE_add) = colnames(DE)
  DE_add[,paste("DEFac", new_group, sep="")] = split_DE[,"DEFacGroup1"]

  updated_X = rbind(X, X_new)
  updated_Y = rbind(Y, Y_new)
  updated_celltype = celltype
  updated_celltype[new_celltype_ind] = new_group
  updated_DE = rbind(DE, DE_add)


  s = Seurat_from_splatter_info(X=updated_X, Y=updated_Y,
                                celltype=updated_celltype,
                                DE=updated_DE)
  s@misc$split = paste(target_group, new_group, sep="+")
  return (s)
}

gene_patterns_from_DE.splatter = function(DE)
{
  binary = 1*(abs(DE - 1) > 1E-8)
  patterns = apply(binary, 1, paste, collapse="")

  return (patterns)
}
