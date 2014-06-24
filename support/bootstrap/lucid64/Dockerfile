FROM ubuntu:lucid

ADD . bootstrap
RUN bootstrap/vagrant/bootstrap-linux.py docker && rm -rf bootstrap
