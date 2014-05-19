#!/bin/sh
# Install buildbot on Ubuntu.

sudo pip -q install buildbot
BUILDBOT_DIR=/var/lib/buildbot
sudo -u vagrant buildbot create-master -r $BUILDBOT_DIR
(crontab -u vagrant -l ||
 echo "@reboot PATH=$PATH:/usr/local/bin buildbot start $BUILDBOT_DIR") |
 crontab -u vagrant -
sudo -u vagrant cp master.cfg $BUILDBOT_DIR/master.cfg
sudo -H -u vagrant buildbot start $BUILDBOT_DIR
