# Clears everything
system(
    "docker container rm rtest-container; docker volume rm rtest-site-library; docker image rm -f rtest-image"
)
# Rebuilds Docker image
system(
    "docker build -t rtest-image . -f tests/docker/Dockerfile"
)

# Install R dependencies
system(
    "docker container kill rtest-container ; docker container rm rtest-container ; docker run --name rtest-container --env R_LIBS=/root/.R/site-library -v rtest-site-library:/root/.R/site-library rtest-image 'cd /package && ./install_deps.sh'"
)

# Run tests
system(
    "docker container kill rtest-container ; docker container rm rtest-container ; docker container run  --rm --name rtest-container  --env R_LIBS=/root/.R/site-library  -v rtest-site-library:/root/.R/site-library rtest-image 'cd /package && ./run_tests.sh'"
)
