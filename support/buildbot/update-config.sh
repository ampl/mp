#!/bin/sh
# This script uploads buildbot config to the server and reconfigs
# the buildbot.

BUILDBOT_DIR=`ssh ampl.com "echo ~buildbot"`
scp master.cfg ampl.com:$BUILDBOT_DIR/master.cfg
ssh -t ampl.com "cd $BUILDBOT_DIR && sudo -H -u buildbot buildbot reconfig"
