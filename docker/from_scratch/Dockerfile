FROM ubuntu:22.04
RUN apt-get update
RUN apt-get install -y software-properties-common dirmngr wget
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get update
ENV TZ=Etc/UTC
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install tzdata
RUN apt-get install -y r-base
RUN apt-get install -y r-base-dev
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
ENV R_HOME=/usr/lib/R
WORKDIR /root
ADD .Rprofile .
ADD .Renviron .

RUN Rscript -e 'install.packages("tidyverse")'