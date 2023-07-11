# Figures - uses Kang celltype datasets
#
# make_signal_figure.manuscript_snr - \sum_i snr_i/1+snr_i for p=0.05, .25, .5 over
# all the cell types over cell specific ISGs

# make_celltype_snr_figure.manuscript_snr - plots signal of top 8 genes for
# the cell types CD14+ Mono, CD8

# make_example_expression_figure.manuscript_snr - shows
# expression in CD14+ Mono for ISG15 and IFIT1

# make_violin_figure.manuscript_snr - shows distribution of snr

make.manuscript_snr = function(m_version)
{
  pdf_f = paste(figure_dir.manuscript(m_version),
                 "/figure_snr_assymetry.pdf", sep="")
  if (!file.exists(pdf_f)) {
    pdf(pdf_f, width=6, height=4)
    g = make_assymetric_signal_figure.manuscript_snr() +
     guides(fill="none", color="none") +
    theme(axis.text=element_text(size=8),
          axis.text.x = element_text(angle = -10, vjust=-.5),
          axis.title=element_text(size=10,face="bold"))
    print(g)
    dev.off()
  }

  pdf_f = paste(figure_dir.manuscript(m_version),
                 "/figure_snr_signal.pdf", sep="")
  if (!file.exists(pdf_f)) {
    pdf(pdf_f, width=6, height=4)
    g = make_signal_figure.manuscript_snr()  +
    theme(axis.text=element_text(size=8),
          axis.text.x = element_text(angle = -10, vjust=-.5),
          axis.title=element_text(size=10,face="bold")) +
    guides(fill="none")
    print(g)
    dev.off()
  }



  pdf_f = paste(figure_dir.manuscript(m_version),
                 "/figure_snr_gene_snr.pdf", sep="")
  if (!file.exists(pdf_f)) {
    pdf(pdf_f, width=6, height=4)
    g = make_gene_snr_figure.manuscript_snr()   +
    theme(axis.text=element_text(size=8),
          axis.text.x = element_text(angle = -10, vjust=-.5),
          axis.title=element_text(size=10,face="bold")) +
    guides(fill="none") + guides(shape="none")

   # theme(axis.text=element_text(size=8), legend.text = element_text(size=8))

    print(g)
    dev.off()
  }

}

make_signal_figure.manuscript_snr = function(min_snr = 0,
                                             overwrite=F)
{
  p_vals = .5

  t = make_table.manuscript_snr() %>%
      dplyr::filter(snr >= min_snr, p==0.5)

  d_plot = plyr::ddply(t, plyr::.(celltype), function(tt) {
       V1 = tt$V1; V2 = tt$V2; mu2 = tt$delta_mu2
       p_signal = sapply(p_vals, function(cp) {
         snr = cp*(1-cp)*mu2/(cp*(1-cp) + cp*V1 + (1-cp)*V2)
         sum(snr/(1+snr))
       })
      data.frame(p=p_vals, signal=p_signal)
  })

  g = ggplot(data=d_plot) +
        geom_bar(aes(x=celltype, y=signal, fill=celltype),
              stat="identity")

  return (g)
}


make_gene_snr_figure.manuscript_snr = function(show_genes=c("CXCL10"),
                                               ngenes=50)
{
  p = 0.5
  celltype = get_celltypes.Kang() %>% setdiff("Kang_all")
  d_all = plyr::adply(celltype, 1, function(cct) {
    print(cct)
    s = permutation.Kang_IFN("none", F, cct, p)
    ISGs = ISG_table.Kang(cct)$gene
    ISGs_active = intersect(ISGs, rownames(s))

    D = s@assays$RNA@data[ISGs_active,]

    d = plyr::adply(1:nrow(D), 1, function(i) {
      g = D[i,]
      gene = rownames(D)[i]

      gs = g[s$stim=="stim"]; gc = g[s$stim=="ctrl"]
      V1 = var(gs); V2 = var(gc)
      dmu2 = (mean(gs) - mean(gc))^2
      data.frame(gene=gene,
                 celltype=cct,
                 p=p,
                 snr=p*(1-p)*dmu2/(p*V1 + (1-p)*V2))
    })

    d = dplyr::arrange(d, desc(snr))
    s_ind = which(is.element(d$gene,show_genes))
    ind = union(1:ngenes, s_ind)
    return (d[ind,])
  })

  d_show = dplyr::filter(d_all, is.element(gene, show_genes))

  g = ggplot() +
    geom_violin(aes(x=celltype, y=snr, fill=celltype), data=d_all) +
    geom_point(aes(x=celltype, y=snr, shape=gene), fill="black",
               data=d_show, size=2)

  return (g)
}


make_assymetric_signal_figure.manuscript_snr = function(min_snr = 0)
{
  p_vals = c(0.05, 0.50, 0.95)

  t = make_table.manuscript_snr() %>%
    dplyr::filter(snr >= min_snr, p==0.5)

  d_plot = plyr::ddply(t, plyr::.(celltype), function(tt) {
    V1 = tt$V1; V2 = tt$V2; mu2 = tt$delta_mu2
    p_signal = sapply(p_vals, function(cp) {
      snr = cp*(1-cp)*mu2/(cp*(1-cp) + cp*V1 + (1-cp)*V2)
      sum(snr/(1+snr))
    })
    data.frame(p=p_vals, signal=p_signal)
  })

  g = ggplot(data=d_plot) +
    geom_bar(aes(x=celltype, y=signal, fill=celltype,
                 color=factor(p)),
             stat="identity", position="dodge")

  return (g)
}



#################################################
make_table.manuscript_snr = function(overwrite=F)
{
  f = paste(dir.Kang_IFN(), "/snr_table.csv", sep="")
  if (!file.exists(f) | overwrite) {
    ct = get_celltypes.Kang() %>% setdiff("Kang_all")

    j = lapply(ct, function(cct) {
      print(cct)
      ISGs = ISG_table.Kang(cct)$gene
      snrs = plyr::adply(c(0.05, 0.25, 0.5), 1, function(p) {

      print(p)
      s = permutation.Kang_IFN("none", F, cct, p)
      cISGs = intersect(ISGs, rownames(s))

      D = s@assays$RNA@data %>% Matrix::as.matrix()
      stim = s$stim == "stim"
      pt = mean(stim)
      delta_mu = (rowMeans(D[cISGs,stim]) - rowMeans(D[cISGs,!stim]))
      delta_mu2 = delta_mu^2

      V1 = apply(D[cISGs,stim], 1, var)
      V2 = apply(D[cISGs,!stim], 1, var)

      snr = pt*(1-pt)*delta_mu2/(pt*V1 + (1-pt)*V2)
      data.frame(celltype=cct, p=p, gene=cISGs, snr=snr,
                   V1=V1, V2=V2, delta_mu2=delta_mu2,
                   delta_mu=delta_mu)
      })
      return (snrs)
    })
    j = do.call(rbind, j)

    write.csv(j, f, row.names = F)
  }
  j = read.csv(f, check.names=F)
  return (j)

}

