#!/bin/bash
Rscript -e 'install.packages("testthat", force=FALSE)'
Rscript -e 'devtools::install_deps(upgrade="always")'