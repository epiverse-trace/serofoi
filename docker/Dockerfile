FROM rocker/tidyverse:4.2.2
VOLUME [ "/root/.R/site-library", "/package" ]
WORKDIR /root
ADD . /package
ADD docker/.bashrc .
# -i ensures that any command will source .bashrc https://stackoverflow.com/a/37286648
ENTRYPOINT ["/bin/bash", "-c", "-i"]