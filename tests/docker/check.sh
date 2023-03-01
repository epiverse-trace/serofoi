#!/bin/bash
R --no-echo --no-restore -e 'devtools::check(manual=FALSE, cran = TRUE, error_on=c("error"))'