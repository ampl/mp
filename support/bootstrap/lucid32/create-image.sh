#!/bin/sh
# Create a base docker image for 32-bit Ubuntu Lucid.

sudo debootstrap --arch=i386 lucid lucid
sudo tar -C lucid -c . | docker import - vitaut/lucid32
