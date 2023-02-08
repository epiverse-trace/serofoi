#!/bin/bash
docker run -v rtest-site-library:/root/.R/site-library --mount src="$(pwd)",target="/package",type=bind -it rtest ./run_tests.sh
