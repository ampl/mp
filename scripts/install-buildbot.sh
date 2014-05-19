#!/bin/sh
# Install buildbot on Ubuntu.

sudo pip install buildbot

BUILDBOT_DIR=/var/lib/buildbot
if id -u buildbot >/dev/null 2>&1; then
  echo buildbot user exists
else
  sudo useradd buildbot
fi
sudo -u buildbot buildbot create-master -r $BUILDBOT_DIR
(crontab -u buildbot -l ||
 echo "@reboot PATH=$PATH:/usr/local/bin buildbot start $BUILDBOT_DIR") |
 crontab -u buildbot -
sudo -u buildbot cp master.cfg $BUILDBOT_DIR/master.cfg
sudo -H -u buildbot buildbot start $BUILDBOT_DIR
