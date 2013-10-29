#!/bin/sh
# Set up buildbot on Mac.

port install buildbot

# Create group buildbot with group id 101.
GROUP=buildbot
dscl . create /groups/$GROUP
dscl . create /groups/$GROUP name $GROUP
dscl . create /groups/$GROUP passwd
dscl . create /groups/$GROUP gid 101

# Create user buildbot.
BUILDBOT_HOME=/var/lib/buildbot
USER=buildbot
mkdir -p $BUILDBOT_HOME
dscl . -create /Users/$USER
dscl . -create /Users/$USER RealName "Buildbot"
dscl . -create /Users/$USER NFSHomeDirectory $BUILDBOT_HOME
dscl . -create /Users/$USER UserShell /bin/bash
dscl . -create /Users/$USER PrimaryGroupID 101
chown $USER:$GROUP $BUILDBOT_HOME

# Hide the user.
defaults write /Library/Preferences/com.apple.loginwindow HiddenUsersList -array-add $USER
