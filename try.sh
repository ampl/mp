#!/bin/bash
# Usage example: ./try.sh -b mp-trusty64 -b mp-win64
ssh -M -S my-ctrl-socket -fnNT -L 4444:localhost:5555 ampl.com
ssh -S my-ctrl-socket -O check ampl.com
buildbot try --connect=pb --master=127.0.0.1:4444 --username=mp --passwd=mp --vc=git $@
ssh -S my-ctrl-socket -O exit ampl.com
