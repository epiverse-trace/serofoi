#!/bin/bash
R --no-echo --no-restore -e 'devtools::check(manual=FALSE, cran = TRUE)'