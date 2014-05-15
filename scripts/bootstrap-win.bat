
cd \vagrant

rem Install Windows SDK.
if exist opt\win64\winsdk if not exist "\Program Files\Microsoft SDKs" (
  opt\win64\winsdk\setup.exe -q -params:ADDLOCAL=ALL
)
