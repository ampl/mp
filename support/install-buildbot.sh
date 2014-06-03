#!/bin/sh
# Install buildbot on Ubuntu.

set -e

sudo pip install buildbot

BUILDBOT_DIR=/var/lib/buildbot
if id -u buildbot >/dev/null 2>&1; then
  echo buildbot user exists
else
  sudo useradd buildbot
fi
sudo mkdir $BUILDBOT_DIR
sudo chown buildbot:buildbot $BUILDBOT_DIR
sudo -u buildbot buildbot create-master -r $BUILDBOT_DIR
(sudo crontab -u buildbot -l ||
 echo "@reboot PATH=$PATH:/usr/local/bin buildbot start $BUILDBOT_DIR") |
 sudo crontab -u buildbot -
sudo -u buildbot cp master.cfg $BUILDBOT_DIR/master.cfg
sudo -H -u buildbot buildbot start $BUILDBOT_DIR
