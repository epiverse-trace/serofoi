# Clear everything
system(
  "docker container rm rtest-container; docker volume rm rtest-site-library; docker image rm -f rtest-image"
)
# Rebuild Docker image
system(
  "docker container kill rtest-container; rm -f inst/extdata/stanmodels/*.rds; docker build -t rtest-image . -f tests/docker/Dockerfile"
)

# Install R dependencies
system(
  "docker container kill rtest-container ; docker container rm rtest-container ; docker run --name rtest-container --env R_LIBS=/root/.R/site-library -v rtest-site-library:/root/.R/site-library rtest-image 'cd /package && ./install_deps.sh'"
)

# Run tests
system(
  "docker container kill rtest-container ; docker container rm rtest-container ; docker container run  --rm --name rtest-container  --env R_LIBS=/root/.R/site-library  -v rtest-site-library:/root/.R/site-library rtest-image 'cd /package && ./run_tests.sh'"
)

# R CMD Check
system(
  "docker container kill rtest-container ; docker container rm rtest-container ; docker container run  --rm --name rtest-container  --env R_LIBS=/root/.R/site-library  -v rtest-site-library:/root/.R/site-library rtest-image 'cd /package && ./check.sh'"
)

# Bash shell
system(
  "docker container run -it  --rm  --env R_LIBS=/root/.R/site-library  -v rtest-site-library:/root/.R/site-library rtest-image bash"
)

# R Shell
system(
  "docker container run -it  --rm  --env R_LIBS=/root/.R/site-library  -v rtest-site-library:/root/.R/site-library rtest-image R"
)
