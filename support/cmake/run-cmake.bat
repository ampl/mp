@echo on
rem Use this script if you want to build with the Windows SDK toolchain.
rem It configures build environment and runs CMake forwarding all command-line
rem arguments to it, for example:
rem
rem   support\cmake\run-cmake -G "Visual Studio 10 Win64" .

set args=%*
set args=%args:"=\"%
for /F "delims=" %%i IN ('cmake "-DARGS=%args%" -P %~dp0/FindSetEnv.cmake') DO (
  call %%i
)
cmake %*
