make_splatter_figure.manuscript_plot = function(m_version)
{
  # full spectrum figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "splatter_spectrum.pdf")
  if (!file.exists(fname)) {
    s = nISG.Kang_IFN(200, "none", T, "splatter_all", 0.5, "joint")
    Y = scaled_counts.dataset(s)
    evals = eigen(t(Y) %*% Y)$values
    d = data.frame(x=evals)

    g = ggplot() + geom_histogram(data=d, aes(x=x),
                                  binwidth=.1,
                                  colour = 1, fill = "white") +
      xlab("eigenvalue") + ylab("count") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))
    pdf(fname, width=3, height=3)
    print(g)
    dev.off()
  }

  # bulk figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "splatter_bulk.pdf")
  if (!file.exists(fname)) {
    s = nISG.Kang_IFN(200, "none", T, "splatter_all", 0.5, "joint")
    s = update_Seurat_splatter_info(s) %>%
      recompute_spike_from_Seurat()

    Xs = s@misc$Xs %>% t
    b = bulk.analysis_util(s, Xs)
    sp = spectrode_MP.analysis_util(b)
    evals = b$vals

    g = plot.spectrode(sp, evals=evals, nbins=50, binwidth=.1) +
        xlim(0,2) +
      xlab("eigenvalue") +
        ylab("count") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=3, height=3)
    print(g)
    dev.off()
  }

  # filtering figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "splatter_filtering.pdf")
  if (!file.exists(fname)) {
    df_f = filtering.manuscript_plot(max_A=3,
                                     permutation_types="none",
                                     datasets="splatter_all")
    g = ggplot() + geom_point(aes(x=A, y=A_spike), size=.8, data=df_f$df_true) +
      geom_line(aes(x=A, y=A_spike), data=df_f$df_pred) +
      xlab("signal strength") + ylab("filtered ss") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=3, height=3)
    print(g)
    dev.off()
  }

  # clustering figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "splatter_clustering.pdf")
  if (!file.exists(fname)) {
    df_c = clustering.manuscript_plot(max_A=6,
                                     permutation_types="none",
                                     datasets="splatter_all")
    g = ggplot() + geom_point(aes(x=A, y=nn), size=.8, data=df_c$df_true) +
      geom_line(aes(x=A, y=nn), data=df_c$df_pred) +
      xlab("signal strength") + ylab("nn metric") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=3, height=3)
    print(g)
    dev.off()
  }

  # umap figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "splatter_umap.pdf")
  if (!file.exists(fname)) {
    d = umap.manuscript_plot("splatter_all", ss=c(0,1,2,3,4,6),
                             permutation_type="none",
                             p=0.5)
    d$label = paste("ss=", round(100*d$A)/100,
                    "  (fss=", round(100*d$A_spike)/100,
                    ",nn=", round(100*d$nn)/100, ")", sep="")
    g = ggplot() + geom_point(aes(x=umap1, y=umap2,
                                   color=cell_state), size=.025, data=d) +
     # geom_text(aes(x=-3.75, y=3, label=round(100*nn)/100), data=d) +
      facet_wrap("label", nrow=3, ncol=2) +
      theme(axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()) + theme_bw()  +
      theme(strip.text.x = element_text(size = 10, face="bold")) +
      guides(color="none")

    pdf(fname, width=5, height=5)
    print(g)
    dev.off()
  }
  return (NULL)
}

make_Kang_celltype.manuscript_plot = function(m_version)
{
  # filtering figure
  datasets = c("FCGR3A+ Monocytes", "CD8 T cells", "CD14+ Monocytes")
  p = c(0.05, 0.25, 0.5)
  df_joint = expand.grid(dataset=datasets, p=p, stringsAsFactors = F)
  datasets = df_joint$dataset
  p = df_joint$p

  # bulk figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "Kang_celltype_bulk.pdf")
  if (!file.exists(fname)) {
    bulk_datasets = setdiff(get_celltypes.Kang(T) %>% dplyr::arrange(n) %>%
                              dplyr::pull(celltype),
                            c("Megakaryocytes", "CD4 T cells", "Kang_all"))
    d = lapply(bulk_datasets, function(dataset) {
      print(dataset)

      s = nISG.Kang_IFN(0, "none", F, dataset, 0.5)
      uv = s@misc$uv
      Xs = uv$u[,1:20] %*% diag(uv$d[1:20]) %*% t(uv$v[,1:20])
      b = bulk.analysis_util(s, Xs)
      sp = spectrode_MP.analysis_util(b)
      # leave out the first because that will be library size correlated
      evals = b$vals[-1]

      param = sp$param
      df_theory = data.frame(x=sp$grid,
                             y=sp$density,
                             dataset=dataset)
      if (param$gamma < 1)
        df_theory$y = df_theory$y/param$gamma
      df_theory$y_scaled = df_theory$y*length(evals)

      df_true = data.frame(x=evals, dataset=dataset)
      return (list(df_theory=df_theory, df_true=df_true))
    })
    df_true = do.call(rbind, lapply(d, function(dd) dd$df_true))
    df_theory = do.call(rbind, lapply(d, function(dd) dd$df_theory))

    df_true$dataset = factor(df_true$dataset, levels=bulk_datasets)
    df_theory$dataset = factor(df_theory$dataset, levels=bulk_datasets)

    # spectrode can't resolve very end
    df_theory = dplyr::filter(df_theory, x > 0.025)

    g = ggplot() +
      geom_histogram(data=df_true, aes(x=x),
                     binwidth=.1,
                     colour = 1, fill = "white") +
      geom_line(data=df_theory, aes(x=x,y=y_scaled*.1),
                color="red", size=1) +
      facet_wrap("dataset", nrow=2, ncol=3) +
      xlim(0,5) +
      xlab("eigenvalue") +
      ylab("count") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=7, height=5)
    print(g)
    dev.off()
  }

  # filtering figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "Kang_celltype_filtering.pdf")
  if (!file.exists(fname)) {

    df_f = filtering.manuscript_plot(max_A=3,
                                     permutation_types="none",
                                     datasets=datasets,
                                     pvals=p)
    df_f$dataset = factor(df_f$dataset, levels=unique(datasets))
    g = ggplot() + geom_point(aes(x=A, y=A_spike), size=.8, data=df_f$df_true) +
      geom_line(aes(x=A, y=A_spike), data=df_f$df_pred) +
      facet_grid(dataset ~ p) +
      xlab("signal strength") + ylab("filtered ss") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=6, height=5)
    print(g)
    dev.off()
  }

  # clustering figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "Kang_celltype_clustering.pdf")
  if (!file.exists(fname)) {
  
    df_c = clustering.manuscript_plot(max_A=6,
                                      permutation_type="none",
                                      datasets=datasets,
                                      pvals=p)
    df_c$dataset = factor(df_c$dataset, levels=unique(datasets))
    g = ggplot() + geom_point(aes(x=A, y=nn), size=.8, data=df_c$df_true) +
      geom_line(aes(x=A, y=nn), data=df_c$df_pred) +
      facet_grid(dataset ~ p) +
      xlab("signal strength") + ylab("nn metric") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=6, height=5)
    print(g)
    dev.off()
  }

  # umap figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "Kang_celltype_umap.pdf")
  if (!file.exists(fname)) {
    d = umap.manuscript_plot(c("CD8 T cells",
                               "FCGR3A+ Monocytes",
                               "CD14+ Monocytes"), ss=c(0,6),
                             permutation_type="none",
                             p=0.5)

    d$dataset = factor(d$dataset, levels=unique(datasets))
    d$label = paste("ss=", d$ss, sep="")

    g = ggplot() + geom_point(aes(x=umap1, y=umap2,
                                  color=cell_state), size=.025, data=d) +
      # geom_text(aes(x=-3.75, y=3, label=round(100*nn)/100), data=d) +
      facet_grid(dataset ~ label) +
      theme(axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()) + theme_bw()  +
      theme(strip.text.x = element_text(size = 10, face="bold")) +
      guides(color="none")

    pdf(fname, width=7, height=7)
    print(g)
    dev.off()
  }
  return (NULL)
}




make_full_dataset.manuscript_plot = function(m_version)
{
  datasets = c(get_all_DE_datasets.dataset_DE(),"Kang_all")
  pt = "none"

  # bulk figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "full_datasets_bulk.pdf")
  if (!file.exists(fname)) {
    d = lapply(datasets, function(dataset) {
      print(dataset)

      if (dataset == "Kang_all")
        s = nISG.Kang_IFN(0, pt, F, dataset, 0.5, "joint")
      else
        s = nISG.dataset_DE(0, pt, dataset)
      uv = s@misc$uv
      Xs = uv$u %*% diag(uv$d) %*% t(uv$v)
      b = bulk.analysis_util(s, Xs)

      B = b$B
      par = param(colSums(B^2), w=1, gamma=ncol(B)/nrow(B))
      sp = spectrode(par, dx=1E-5, epsilon=1E-8)

      # leave out the first because that will be library size correlated
      evals = b$vals[-1]

      param = sp$param
      df_theory = data.frame(x=sp$grid,
                             y=sp$density,
                             dataset=dataset)
      if (param$gamma < 1)
        df_theory$y = df_theory$y/param$gamma
      df_theory$y_scaled = df_theory$y*length(evals)

      df_true = data.frame(x=evals, dataset=dataset)
      return (list(df_theory=df_theory, df_true=df_true))
    })
    df_true = do.call(rbind, lapply(d, function(dd) dd$df_true))
    df_theory = do.call(rbind, lapply(d, function(dd) dd$df_theory))

    df_true$dataset = fancy_names.manuscript_plot(df_true$dataset)
    df_theory$dataset = fancy_names.manuscript_plot(df_theory$dataset)

    # spectrode can't resolve very end
    df_theory = dplyr::filter(df_theory, x > 0.025)

    g = ggplot() +
      geom_histogram(data=df_true, aes(x=x),
                     binwidth=.1,
                     colour = 1, fill = "white") +
      geom_line(data=df_theory, aes(x=x,y=y_scaled*.1),
                color="red", size=1) +
      facet_wrap("dataset", nrow=2, ncol=2) +
      xlim(0,2.25) +
      xlab("eigenvalue") +
      ylab("count") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=5, height=5)
    print(g)
    dev.off()
  }


  # filtering figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "full_datasets_filtering.pdf")
  if (!file.exists(fname)) {

    df_f = filtering.manuscript_plot(max_A=3,
                                     permutation_types=pt,
                                     datasets=datasets)
    df_f$df_true$dataset = fancy_names.manuscript_plot(df_f$df_true$dataset)
    df_f$df_pred$dataset = fancy_names.manuscript_plot(df_f$df_pred$dataset)
    g = ggplot() + geom_point(aes(x=A, y=A_spike), size=.8, data=df_f$df_true) +
      geom_line(aes(x=A, y=A_spike), data=df_f$df_pred) +
      facet_wrap("dataset", nrow=2, ncol=2) +
      xlab("signal strength") + ylab("filtered ss") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=6, height=6)
    print(g)
    dev.off()
  }

  # clustering figure
  fname = here::here(figure_dir.manuscript(m_version),
                     "full_datasets_clustering.pdf")
  if (!file.exists(fname)) {
    df_c = clustering.manuscript_plot(max_A=8.5,
                                      permutation_type=pt,
                                      datasets=datasets)
    df_I = clustering.manuscript_plot(max_A=8.5, "module", "Intestine")

    df_c$df_true$dataset = fancy_names.manuscript_plot(df_c$df_true$dataset)
    df_c$df_pred$dataset = fancy_names.manuscript_plot(df_c$df_pred$dataset)
    df_I$df_true$dataset = fancy_names.manuscript_plot(df_I$df_true$dataset)

    g = ggplot() + geom_point(aes(x=A, y=nn), size=.8, data=df_c$df_true) +
      geom_line(aes(x=A, y=nn), data=df_c$df_pred) +
      geom_point(aes(x=A, y=nn), size=.8, col="blue", shape="triangle", data=df_I$df_true) +
      facet_wrap("dataset", ncol=2, nrow=2) +
      xlab("signal strength") + ylab("nn metric") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"))

    pdf(fname, width=6, height=6)
    print(g)
    dev.off()
  }


  # umap figure
  fname_stim = here::here(figure_dir.manuscript(m_version),
                     "full_datasets_umap_stim.pdf")
  #fname_celltype = here::here(figure_dir.manuscript(m_version),
  #                        "full_datasets_umap_celltype.pdf")
  if (!file.exists(fname_stim)) {
    d1 = umap.manuscript_plot(c("Intestine"), ss=c(0,5),
                             permutation_type="none")
    d2 = umap.manuscript_plot(c("Intestine"), ss=c(0,5),
                              permutation_type="module")
    d1$dataset = fancy_names.manuscript_plot(d1$dataset)
    d2$dataset = "permuted"

    d = rbind(d1, d2)
    d$label = paste("ss=",
                    ifelse(d$A < 2, 0, 5), sep="")

    g_stim = ggplot() + geom_point(aes(x=umap1, y=umap2,
                                  color=cell_state), size=.025, data=d) +
      # geom_text(aes(x=-3.75, y=3, label=round(100*nn)/100), data=d) +
      facet_grid(dataset ~ label) +
      theme(axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()) + theme_bw()  +
      theme(strip.text.x = element_text(size = 10, face="bold")) +
      guides(color="none")

    # g_celltype = ggplot() + geom_point(aes(x=umap1, y=umap2,
    #                                    color=celltype), size=.025, data=d) +
    #   # geom_text(aes(x=-3.75, y=3, label=round(100*nn)/100), data=d) +
    #   facet_wrap("label", nrow=1) +
    #   theme(axis.title=element_blank(),
    #         axis.text=element_blank(),
    #         axis.ticks=element_blank()) + theme_bw()  +
    #   theme(strip.text.x = element_text(size = 10, face="bold")) +
    #   guides(color = guide_legend(override.aes = list(size=5)))


    pdf(fname_stim, width=6, height=6)
    print(g_stim)
    dev.off()
    #pdf(fname_celltype, width=6, height=3)
    #print(g_celltype)
    #dev.off()
  }
  return (NULL)
}



############################################################
fancy_names.manuscript_plot = function(datasets)
{
  datasets = as.character(datasets)
  convert = sapply(datasets, function(d) {
    switch(d,
           Kang_all="Kang et al.",
           Intestine="Mead et al.",
           HumanNasal="Ordovas-Montanes et al.",
           Zheng="Zheng et al.",
           stop(problem))
  })

  return (convert)
}

filtering.manuscript_plot = function(max_A=4,
                                     permutation_types,
                                     datasets,
                                     pvals=NULL)
{
  if (is.null(pvals))
    pvals = rep(0.5, length(datasets))
  df_true = plyr::adply(1:length(datasets), 1, function(i) {
    dd = datasets[i]
    if (dd %in% get_all_DE_datasets.dataset_DE())
      df = filtering_table.dataset_DE(dd)
    else
      df = filtering_table.manuscript_filtering() %>%
        dplyr::filter(celltype==dd, p==pvals[i]) %>% dplyr::mutate(dataset=celltype)

    if (dd %in% get_all_DE_datasets.dataset_DE())
      df = dplyr::select(df, dataset, A, A_spike,
                         permutation_type)
    else
      df = dplyr::select(df, dataset, A, A_spike,
                         permutation_type, p)
  })
  df_pred = plyr::adply(1:length(datasets), 1, function(i) {
    dd = datasets[i]
    if (dd %in% get_all_DE_datasets.dataset_DE())
      df = filtering_theory_table.dataset_DE(dd)
    else {
      df = filtering_theory_table.manuscript_filtering() %>%
        dplyr::filter(celltype==dd, p==pvals[i]) %>% dplyr::mutate(dataset=celltype)
    }

    if (dd %in% get_all_DE_datasets.dataset_DE())
      df = dplyr::select(df, dataset, A, A_spike,
                         permutation_type)
    else
      df = dplyr::select(df, dataset, A, A_spike,
                         permutation_type, p)
  })

  udatasets = unique(datasets)
  df_true$dataset = factor(df_true$dataset, levels = udatasets)
  df_pred$dataset = factor(df_pred$dataset, levels = udatasets)

  df_true = dplyr::filter(df_true, A <= max_A)
  df_pred = dplyr::filter(df_pred, A <= max_A)

  df_true = df_true[is.element(df_true$permutation_type, permutation_types),]
  df_pred = df_pred[is.element(df_pred$permutation_type, permutation_types),]

  return (list(df_true=df_true, df_pred=df_pred))
}


clustering.manuscript_plot = function(max_A=4, permutation_types,
                                      datasets, pvals=NULL)
{
  if (is.null(pvals))
    pvals = rep(0.5, length(datasets))

  df_true = plyr::adply(1:length(datasets), 1, function(i) {
    dd = datasets[i]
    print(dd)
    if (dd %in% get_all_DE_datasets.dataset_DE())
      df = clustering_table.manuscript_DE(dd)
    else
      df = clustering_table.manuscript_clustering() %>%
        dplyr::filter(celltype==dd, p==pvals[i]) %>% dplyr::mutate(dataset=celltype)

    if (dd %in% get_all_DE_datasets.dataset_DE())
      df = dplyr::select(df, dataset, A, nn,
                         permutation_type)
    else
      df = dplyr::select(df, dataset, A, nn,
                         permutation_type, p)
  })
  df_pred = plyr::adply(1:length(datasets), 1, function(i) {
    dd = datasets[i]
    print(dd)
    if (dd %in% get_all_DE_datasets.dataset_DE())
      df = clustering_table_theory.manuscript_DE(dd)
    else {
      df = plyr::adply(c( "none"), 1, function(pt)
        clustering_table_theory.manuscript_clustering(pt)) %>%
        dplyr::filter(celltype==dd, p==pvals[i]) %>% dplyr::mutate(dataset=celltype)
    }

    if (dd %in% get_all_DE_datasets.dataset_DE())
      df = dplyr::select(df, dataset, A, nn,
                         permutation_type)
    else
      df = dplyr::select(df, dataset, A, nn,
                         permutation_type, p)
  })

  udatasets = unique(datasets)
  df_true$dataset = factor(df_true$dataset, levels = udatasets)
  df_pred$dataset = factor(df_pred$dataset, levels = udatasets)

  df_true = dplyr::filter(df_true, A <= max_A)
  df_pred = dplyr::filter(df_pred, A <= max_A)


  df_true = df_true[is.element(df_true$permutation_type, permutation_types),]
  df_pred = df_pred[is.element(df_pred$permutation_type, permutation_types),]

  return (list(df_true=df_true, df_pred=df_pred))
}

umap.manuscript_plot = function(datasets, ss,
                permutation_type="none",
                p=0.5)
{
  d_umap = plyr::adply(datasets, 1, function(dataset) {
    if (dataset %in% get_all_DE_datasets.dataset_DE()) {
      dc = clustering_table.manuscript_DE(dataset)
      df = filtering_table.dataset_DE(dataset)
    }
    else {
      dc = clustering_table.manuscript_clustering() %>%
        dplyr::filter(celltype==dataset, p==!!p) %>% dplyr::mutate(dataset=celltype)
      df = filtering_table.manuscript_filtering() %>%
        dplyr::filter(celltype==dataset, p==!!p) %>% dplyr::mutate(dataset=celltype)
    }

    dc = dplyr::filter(dc, permutation_type==!!permutation_type)
    df = dplyr::filter(df, permutation_type==!!permutation_type)
    if (!all(dc$A == df$A))
      stop("problem")
    d = dplyr::mutate(df, nn=dc$nn)

    nISG = sapply(ss, function(css) {
      ind = which.min(abs(d$A - css))
      print(d[ind,])
      return (d$nISG[ind])
    })
    dd_umap = plyr::adply(1:length(nISG), 1, function(i) {

      if (dataset %in% get_all_DE_datasets.dataset_DE())
        s = nISG.dataset_DE(nISG[i], permutation_type, dataset)
      else {
        if (grepl("all", dataset))
          ISG_type = "joint"
        else
          ISG_type = NULL
        if (dataset %in% get_celltypes.Kang())
          insilico = F
        else
          insilico = T
        s = nISG.Kang_IFN(nISG[i], permutation_type, insilico,
                          dataset, p, ISG_type)
      }

      dn = d[d$nISG==nISG[i],]
      nn =  dn$nn
      A_spike = dn$A_spike

      snr = s@misc$snr
      df = data.frame(umap1=s@reductions$umap@cell.embeddings[,1],
                    umap2=s@reductions$umap@cell.embeddings[,2],
                    dataset=dataset,
                    nn =nn,
                    A_spike=A_spike,
                    A = sum(snr/(1+snr)),
                    ss=ss[i],
                    index=i,
                    cell_state=s$stim,
                    celltype=s$celltype)

      return (df)
    })
  })

  return (d_umap)
}
