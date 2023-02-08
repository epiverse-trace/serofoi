#!/bin/bash
Rscript -e 'install.packages("testthat")'
Rscript -e 'install.packages("covr")'
# Rscript -e 'devtools::install(upgrade=TRUE)'
Rscript -e 'devtools::install_deps(upgrade="always")'