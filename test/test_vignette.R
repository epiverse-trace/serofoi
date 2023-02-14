# library(remotes)
# remotes::install_github("TRACE-LAC/serofoi", ref = "dev")
library(serofoi)


library(usethis)
# usethis::use_pkgdown()
# pkgdown::build_site()
use_github_pages(branch = "dev-zulma-vignette", path = "/", cname = NA)
