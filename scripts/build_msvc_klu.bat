@echo off
setlocal

set "DFN_BUILD_DIR=build-klu"
set "VCPKG_TOOLCHAIN=%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake"

if not exist "%VCPKG_TOOLCHAIN%" (
  echo error: vcpkg toolchain not found under %VCPKG_ROOT%
  echo Set VCPKG_ROOT to the vcpkg checkout.
  exit /b 1
)

call "%~dp0build_msvc.bat" ^
  -DCMAKE_TOOLCHAIN_FILE=%VCPKG_TOOLCHAIN% ^
  -DVCPKG_TARGET_TRIPLET=x64-windows ^
  -DVCPKG_OVERLAY_PORTS=%~dp0..\cmake\vcpkg-overlay ^
  -DDFN_USE_KLU=ON
