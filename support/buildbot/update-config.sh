#!/bin/bash
# This script uploads buildbot config to the server and reconfigs
# the buildbot.
BASEDIR=`dirname $0`
cd $BASEDIR

if [ $# -ne 1 ]; then
    echo -e "Usage:"
    echo -e "  $0 <password>"
    exit 1
else
    PASSWORD=$1
fi

BUILDBOT_DIR=`ssh ampl.com "echo ~buildbot"`
cp master.cfg tmp.cfg
sed -r 's/\[\("admin","password"\)\]/\[\("admin","'$PASSWORD'"\)\]/g' \
    master.cfg > tmp.cfg
scp tmp.cfg ampl.com:$BUILDBOT_DIR/master.cfg
scp templates/layout.html ampl.com:$BUILDBOT_DIR/templates/layout.html
ssh -t ampl.com "cd $BUILDBOT_DIR; sudo -H -u buildbot buildbot reconfig"
rm tmp.cfg
