
cd \vagrant

rem Install .NET Framework 4 for msbuild.
rem This requires vagrant-windows plugin version 1.7.0.pre.2 or later.
rem See https://github.com/WinRb/vagrant-windows/pull/189
if exist opt\win64\dotNetFx40_Full_x86_x64.exe (
  opt\win64\dotNetFx40_Full_x86_x64.exe /q
)

rem Install Windows SDK.
if exist opt\win64\winsdk if not exist "\Program Files\Microsoft SDKs" (
  opt\win64\winsdk\setup.exe -q -params:ADDLOCAL=ALL
)

rem Install 64-bit JDK.
if exist opt\win64\jdk-7u55-windows-x64.exe (
  opt\win64\jdk-7u55-windows-x64.exe /s
)

rem Install 32-bit JDK.
if exist opt\win64\jdk-7u55-windows-i586.exe (
  opt\win64\jdk-7u55-windows-i586.exe /s
)

rem Install CMake.
bitsadmin /transfer cmake %
  http://www.cmake.org/files/v2.8/cmake-2.8.12.2-win32-x86.zip %
  \vagrant\cmake.zip

rem TODO: install CMake, mingw
