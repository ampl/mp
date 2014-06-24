#!/bin/sh
# Create a base docker image for 32-bit Ubuntu Lucid.

set -e
sudo debootstrap --arch=i386 lucid lucid
sudo tar -C lucid -c . | docker import - lucid32
docker run -v `pwd`/..:/support --name lucid32 lucid32 \
  support/bootstrap/bootstrap-linux.py docker
docker commit lucid32 vitaut/ampl:lucid32
