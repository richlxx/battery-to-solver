@echo off
setlocal

set "DFN_BUILD_DIR=build-klu"
set "VCPKG_INSTALLED=%VCPKG_ROOT%\installed\x64-windows"

if not exist "%VCPKG_INSTALLED%\share\sundials\SUNDIALSConfig.cmake" (
  echo error: SUNDIALS with KLU is not installed under %VCPKG_INSTALLED%
  echo Run: vcpkg install --classic --overlay-ports=cmake/vcpkg-overlay "sundials[klu]:x64-windows"
  exit /b 1
)

call "%~dp0build_msvc.bat" -DDFN_USE_KLU=ON -DCMAKE_PREFIX_PATH=%VCPKG_INSTALLED%
