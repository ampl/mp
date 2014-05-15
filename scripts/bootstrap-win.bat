
cd \vagrant

rem Install .NET Framework 4 for msbuild.
rem This requires vagrant-windows plugin version 1.7.0.pre.2 or later.
rem See https://github.com/WinRb/vagrant-windows/pull/189
if not exist "\Windows\Microsoft.NET\Framework64\v4.0.30319" (
  if exist opt\win64\dotNetFx40_Full_x86_x64.exe (
    opt\win64\dotNetFx40_Full_x86_x64.exe /q /norestart
  )
)

rem Install Windows SDK.
if not exist "\Program Files\Microsoft SDKs" if exist opt\win64\winsdk (
  opt\win64\winsdk\setup.exe -q -params:ADDLOCAL=ALL
)

rem Install 64-bit JDK.
if not exist "\Program Files\Java\jdk1.7.0_55" (
  if exist opt\win64\jdk-7u55-windows-x64.exe (
    opt\win64\jdk-7u55-windows-x64.exe /s
  )
)

rem Install 32-bit JDK.
if not exist "\Program Files (x86)\Java\jdk1.7.0_55" (
  if exist opt\win32\jdk-7u55-windows-i586.exe (
    opt\win32\jdk-7u55-windows-i586.exe /s
  )
)

rem Install CMake.
if not exist "\Program Files\CMake" (
  bitsadmin /transfer cmake ^
    http://www.cmake.org/files/v2.8/cmake-2.8.12.2-win32-x86.zip ^
    C:\vagrant\cmake.zip
  "\Program Files (x86)\Java\jdk1.7.0_55\bin\jar" xf cmake.zip
  move cmake-2.8.12.2-win32-x86 "\Program Files\CMake"
  setx Path "%Path%;C:\Program Files\CMake\bin"

  if exist opt\win64\root (
    copy opt\win64\root \
  )
)
