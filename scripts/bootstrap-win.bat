
cd \vagrant

rem Install .NET Framework 4 for msbuild.
if exist opt\win64\dotNetFx40_Full_x86_x64.exe (
  opt\win64\dotNetFx40_Full_x86_x64.exe /q
)

rem Install Windows SDK.
if exist opt\win64\winsdk if not exist "\Program Files\Microsoft SDKs" (
  opt\win64\winsdk\setup.exe -q -params:ADDLOCAL=ALL
)

rem TODO: install Java, CMake, mingw

rem Install CMake.
bitsadmin /transfer cmake %
  http://www.cmake.org/files/v2.8/cmake-2.8.12.2-win32-x86.zip %
  \vagrant\cmake.zip
