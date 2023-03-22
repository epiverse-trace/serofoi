FROM rocker/tidyverse:4.2.2
VOLUME [ "/root/.R/site-library"]
ARG docker_folder="tests/docker"
WORKDIR /package
ADD . .
ADD ${docker_folder}/*.sh ./
# ADD ${docker_folder}/install_deps.sh .

WORKDIR /root
ADD ${docker_folder}/.bashrc .
# -i ensures that any command will source .bashrc https://stackoverflow.com/a/37286648
# ENTRYPOINT ["/bin/bash", "-c", "-i"]
ENTRYPOINT ["/bin/bash", "-c"]