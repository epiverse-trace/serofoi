# TODO Convert to testthat. Move to `tests/`
library(remotes)
remotes::install_github("TRACE-LAC/serofoi",
                        ref = "dev",
                        force = TRUE)


library(serofoi)
library(usethis)
library(roxygen2)

devtools::document()
usethis::use_pkgdown_github_pages()
pkgdown::build_site()
use_github_pages(branch = "dev-webrd", path = "/", cname = NA)

