# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.define "lucid32" do |lucid32|
    lucid32.vm.box = "lucid32"
    lucid32.vm.box_url = "http://files.vagrantup.com/lucid32.box"
    lucid32.vm.synced_folder "scripts/vagrant/lucid32/archives", "/var/cache/apt/archives"
    lucid32.vm.provision :shell, :path => "scripts/bootstrap-ubuntu.sh"
  end

  config.vm.define "lucid64", primary: true do |lucid64|
    lucid64.vm.box = "lucid64"
    lucid64.vm.box_url = "http://files.vagrantup.com/lucid64.box"
    lucid64.vm.synced_folder "scripts/vagrant/lucid64/archives", "/var/cache/apt/archives"
    lucid64.vm.provision :shell, :path => "scripts/bootstrap-ubuntu.sh"
  end

  config.vm.define "macosx" do |macosx|
    macosx.vm.box = "macosx"
  end
end
