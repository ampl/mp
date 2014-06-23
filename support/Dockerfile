FROM ubuntu:lucid

ADD . bootstrap
RUN bootstrap/vagrant/bootstrap-linux.py buildbot && rm -rf bootstrap
