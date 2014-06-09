rem Download Python to get sane command-line environment
rem and run the main bootstrap script.

rem If we are in a VM managed by Vagrant, then do everything in
rem the shared /vagrant directory to avoid growth of the VM drive.
if exist \vagrant cd \vagrant

cscript /nologo scripts\wget.js ^
  http://download.microsoft.com/download/d/2/4/d242c3fb-da5a-4542-ad66-f9661d0a8d19/vcredist_x64.exe ^
  vcredist_x64.exe
vcredist_x64.exe /q:a
del vcredist_x64.exe
cscript /nologo scripts\wget.js ^
  https://www.python.org/ftp/python/2.7.6/python-2.7.6.amd64.msi python.msi
msiexec /i python.msi ALLUSERS=1
\Python27\python support\vagrant\bootstrap-windows.py
