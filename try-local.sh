#!/bin/bash
# Usage example: ./try-local.sh -b mp-trusty64 -b mp-win64
buildbot try --connect=pb --master=127.0.0.1:5555 --username=mp --passwd=mp --vc=git $@

