# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.box = "lucid32"
  config.vm.synced_folder "scripts/vagrant/lucid32/archives", "/var/cache/apt/archives"
#  config.vm.provision :shell, :path => "scripts/bootstrap-ubuntu.sh"
end
