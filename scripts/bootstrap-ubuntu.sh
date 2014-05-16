#!/bin/sh
# Set up build environment on Ubuntu.

set -e

CMAKE=cmake-2.8.12.2-Linux-i386
if [ -e /opt/$CMAKE ]; then exit; fi

# If we are in a VM managed by Vagrant, then do everything in the shared
# /vagrant directory to avoid growth of the VM drive.
if [ -e /vagrant ]; then
  cd /vagrant
  # Copy optional dependencies.
  if [ -e opt ]; then
    sudo cp -r opt/linux-`uname -m`/* /opt
  fi
fi

sudo apt-get update

# Install build tools.
sudo apt-get install -y git-core gcc g++ gfortran ccache \
     make python-pip python-dev default-jdk unixodbc-dev

# Install CMake.
wget -nv http://www.cmake.org/files/v2.8/$CMAKE.tar.gz
tar xzf $CMAKE.tar.gz
sudo rm -rf $CMAKE.tar.gz /opt/$CMAKE
sudo mv $CMAKE /opt
sudo ln -sf /opt/$CMAKE/bin/cmake /usr/local/bin

# Set up ccache links.
sudo ln -sf /usr/bin/ccache /usr/local/bin/gcc
sudo ln -sf /usr/bin/ccache /usr/local/bin/cc
sudo ln -sf /usr/bin/ccache /usr/local/bin/g++
sudo ln -sf /usr/bin/ccache /usr/local/bin/c++

# Install f90cache.
rm -rf f90cache
git clone https://github.com/vitaut/f90cache.git
cd f90cache
./configure
make
sudo make install
cd ..
rm -rf f90cache
sudo ln -sf /usr/local/bin/f90cache /usr/local/bin/gfortran-4.4

# Install buildbot.
if [ `uname -m` = "x86_64" ]; then
  sudo pip -q install buildbot
  BUILDBOT_BASEDIR=/home/vagrant/master
  sudo -u vagrant buildbot create-master -r $BUILDBOT_BASEDIR
  (crontab -u vagrant -l ||
   echo "@reboot PATH=$PATH:/usr/local/bin buildbot start $BUILDBOT_BASEDIR") |
   crontab -u vagrant -
  sudo -u vagrant cp scripts/master.cfg $BUILDBOT_BASEDIR/master.cfg
  sudo -u vagrant buildbot start $BUILDBOT_BASEDIR
fi
sudo pip -q install buildbot-slave
