#!/bin/bash
docker build -t rtest . --no-cache  --build-arg MAKEFLAGS="${MAKEFLAGS}" -f docker/Dockerfile