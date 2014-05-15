
cd /vagrant

rem Install Windows SDK.
if exist opt/win64/winsdk (
  opt/win64/winsdk/setup.exe -q -params:ADDLOCAL=ALL
)
