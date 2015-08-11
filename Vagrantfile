# -*- mode: ruby -*-
# vi: set ft=ruby :
#
# Vagrant configuration for AMPL/MP.
#
# Usage:
#   > vagrant up lucid64
#   ...
#   lucid64: Container created: <container-id>
#   ...
#   > docker exec -it <container-id> bash

require 'pathname'

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

# Path to directory containing optional dependencies.
OPT_DIR = ENV["AMPL_OPT_DIR"]

def configure_docker(config, image, arch)
  config.vm.provider "docker" do |d|
    d.image = "vitaut/ampl:" + image
    if OPT_DIR
      dir = OPT_DIR + "/linux-" + arch + "/*/"
      d.volumes = Pathname.glob(dir).map { |p| p.to_s + ":/opt/" + p.basename.to_s }
    end
    d.cmd = ["sudo", "-H", "-u", "buildbot", "buildslave", "start",
             "--nodaemon", "/var/lib/buildbot/slave"]
  end
end

def configure_virtualbox(config, cpus, memory, port)
  # This requires VirtualBox Extension Pack to be installed on the host.
  config.vm.provider "virtualbox" do |v|
    v.cpus = cpus
    v.memory = memory
    v.customize ["modifyvm", :id, "--vrde", "on", "--vrdeauthtype", "external",
                                  "--vrdeport", port.to_s]
  end
end
  
Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  # Linux boxes don't use provisioning. To update them, use
  # support/bootstrap/create-lucid-image.py script.
  config.vm.define "lucid32" do |c|
    configure_docker(c, "lucid32", "i686")
  end

  config.vm.define "lucid64", primary: true do |c|
    configure_docker(c, "lucid64", "x86_64")
  end

  config.vm.define "osx-ml" do |c|
    configure_virtualbox(config, 2, 1024, 5000)
    c.vm.box = "osx-ml"
    c.vm.network :private_network, ip: "10.11.12.13"
    # Options "nolocks" and "locallocks" are required for mounting DMG files
    # from an NFS share.
    mount_options = ["resvport", "nolocks", "locallocks"]
    c.vm.synced_folder ".", "/vagrant", :type => "nfs", :mount_options => mount_options
    if OPT_DIR
      c.vm.synced_folder OPT_DIR, "/mnt/opt", :type => "nfs", :mount_options => mount_options
    end
    c.vm.provision :shell, :inline => "/vagrant/support/bootstrap/bootstrap-osx.py"
  end

  config.vm.define "win2008" do |c|
    configure_virtualbox(config, 4, 2048, 5000)
    c.vm.box = "win2008"
    c.vm.guest = :windows
    c.vm.communicator = "winrm"
    # Write the output to provision.log because of the issue
    # https://github.com/mitchellh/vagrant/issues/3866
    c.vm.provision "shell",
      inline: "\\vagrant\\support\\bootstrap\\bootstrap-windows.bat " +
              "> \\vagrant\\support\\bootstrap\\provision.log 2>&1"
  end
end
