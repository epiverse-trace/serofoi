#!/bin/bash
docker run -v rtest-site-library:/root/.R/site-library --mount src="$(pwd)",target="/package" -it rtest