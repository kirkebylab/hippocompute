suppressMessages(library(rsconnect))
# deploy/deploy-shinyapps.R
# usethis::use_build_ignore("deploy")
rsconnect::setAccountInfo(
  Sys.getenv("SHINYAPPS_ACCOUNT"), # repo secret in GitHub settings
  Sys.getenv("SHINYAPPS_TOKEN"),   # repo secret in GitHub settings
  Sys.getenv("SHINYAPPS_SECRET")   # repo secret in GitHub settings
)
rsconnect::deployApp(
  appName = "hippocompute",
  # exclude hidden files and renv directory (if present)
  appFiles = setdiff(list.files(), "renv")
)
