figure_dir.manuscript = function(submission_version=NULL)
{
  f = here::here(paste("figures", submission_version, sep=""))
  if (!dir.exists(f))
    dir.create(f)

  return (f)
}
