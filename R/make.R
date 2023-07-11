make_everything = function()
{
  make_all_datasets()
  make_perturbation_datasets()
  make_analysis()
  make_figures()
}

make_all_datasets = function()
{
  # create the datasets
  unprocessed.Zheng()
  unprocessed.HumanNasal()
  unprocessed.Intestine()
  make_all.Kang()
  make.Kang_insilico()
}

make_perturbation_datasets = function()
{
  analysis_dir = here::here("analysis")
  if (!dir.exists(analysis_dir))
    dir.create(analysis_dir)
  make.dataset_DE()
  make.Kang_IFN()
}

make_analysis = function()
{
  # make the non-Kang et al. dataset analysis
  make.manuscript_DE()
  # make the Kang et al. analysis
  bulk_profile_table.analysis_util()
  make.manuscript_filtering()
  make.manuscript_clustering()
}

make_figures = function()
{
  make_splatter_figure.manuscript_plot(NULL)
  make_Kang_celltype.manuscript_plot(NULL)
  make_full_dataset.manuscript_plot(NULL)
  make.manuscript_snr(NULL)
  make.manuscript_odds_and_ends(NULL)
}