#!/bin/bash
# Script to update the virtual machines that are used as buildslaves

cd `dirname $0`
CUR_MIRROR="http://archive.ubuntu.com/ubuntu"
OLD_MIRROR="http://old-releases.ubuntu.com/ubuntu"
BOOTSTRAP_LINUX="support/bootstrap/create-ubuntu-image.py"

while true;
do
    case "$1" in
    osx-ml)
        vagrant destroy $1
        vagrant up $1
        shift 1;;

    win2008)
        vagrant destroy $1
        vagrant up $1        
        shift 1;;     

    lucid32)
        vagrant destroy $1
        sudo docker rmi buildbot:$1
        sudo python $BOOTSTRAP_LINUX 32 lucid $1 $OLD_MIRROR || exit 1 
        vagrant up $1
        shift 1;;
    lucid64)
        vagrant destroy $1
        sudo docker rmi buildbot:$1
        sudo python $BOOTSTRAP_LINUX 64 lucid $1 $OLD_MIRROR || exit 1 
        vagrant up $1
        shift 1;;

    precise32)
        vagrant destroy $1
        sudo docker rmi buildbot:$1
        sudo python $BOOTSTRAP_LINUX 32 precise $1 $CUR_MIRROR || exit 1 
        vagrant up $1
        shift 1;;
    precise64)
        vagrant destroy $1
        sudo docker rmi buildbot:$1
        sudo python $BOOTSTRAP_LINUX 64 precise $1 $CUR_MIRROR || exit 1 
        vagrant up $1
        shift 1;;

    trusty32)
        vagrant destroy $1
        sudo docker rmi buildbot:$1
        sudo python $BOOTSTRAP_LINUX 32 trusty $1 $CUR_MIRROR || exit 1 
        vagrant up $1
        shift 1;;
    trusty64)
        vagrant destroy $1
        sudo docker rmi buildbot:$1
        sudo python $BOOTSTRAP_LINUX 64 trusty $1 $CUR_MIRROR || exit 1 
        vagrant up $1
        shift 1;;

    xenial32)
        vagrant destroy $1
        sudo docker rmi buildbot:$1
        sudo python $BOOTSTRAP_LINUX 32 xenial $1 $CUR_MIRROR || exit 1 
        vagrant up $1
        shift 1;;
    xenial64)
        vagrant destroy $1
        sudo docker rmi buildbot:$1
        sudo python $BOOTSTRAP_LINUX 64 xenial $1 $CUR_MIRROR || exit 1 
        vagrant up $1
        shift 1;;

    proxy)
        docker stop squid
        docker rm squid
        docker pull sameersbn/squid:3.3.8-19
        docker run --name squid -d \
            --publish 3128:3128 \
            --volume `pwd`/support/squid.conf:/etc/squid3/squid.conf \
            --volume /srv/docker/squid/cache:/var/spool/squid3 \
            sameersbn/squid:3.3.8-19
        shift 1;;

    *)
        if [[ -n "$1" ]]; then
            echo "Invalid option!"
            exit 1
        else
            break
        fi        
    esac
done
