FROM rocker/tidyverse:4.2.2
WORKDIR /serofoi
ADD . .
RUN Rscript -e 'devtools::install(upgrade=TRUE)'
CMD [ "/serofoi/run_tests.sh" ]