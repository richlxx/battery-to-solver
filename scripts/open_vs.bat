@echo off
setlocal

cd /d "%~dp0.."

set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" (
  echo error: vswhere.exe not found
  exit /b 1
)

for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSPATH=%%i"
if not defined VSPATH (
  echo error: Visual Studio with C++ tools not found
  exit /b 1
)

set "CMAKE_EXE=%VSPATH%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if not exist "%CMAKE_EXE%" (
  echo error: CMake bundled with Visual Studio not found
  exit /b 1
)

if not defined DFN_VS_DIR set "DFN_VS_DIR=build-vs"

if exist "%DFN_VS_DIR%\CMakeCache.txt" goto :build

rem First configure: enable the SUNDIALS/KLU backend when it is already
rem installed, so that dfn-sim (the startup project) can actually run.
set "SUNDIALS_PREFIX="
if defined VCPKG_ROOT if exist "%VCPKG_ROOT%\installed\x64-windows\share\sundials\SUNDIALSConfig.cmake" set "SUNDIALS_PREFIX=%VCPKG_ROOT%\installed\x64-windows"
if not defined SUNDIALS_PREFIX if exist "%CD%\vcpkg_installed\x64-windows\share\sundials\SUNDIALSConfig.cmake" set "SUNDIALS_PREFIX=%CD%\vcpkg_installed\x64-windows"
if not defined SUNDIALS_PREFIX if exist "%CD%\build-klu\vcpkg_installed\x64-windows\share\sundials\SUNDIALSConfig.cmake" set "SUNDIALS_PREFIX=%CD%\build-klu\vcpkg_installed\x64-windows"

if defined SUNDIALS_PREFIX (
  echo Configuring with SUNDIALS/KLU: %SUNDIALS_PREFIX%
  "%CMAKE_EXE%" -S . -B "%DFN_VS_DIR%" -G "Visual Studio 17 2022" -A x64 ^
    -DDFN_USE_KLU=ON -DCMAKE_PREFIX_PATH="%SUNDIALS_PREFIX%"
  goto :check
)

if defined VCPKG_ROOT (
  echo SUNDIALS is not installed yet; letting vcpkg install the manifest.
  "%CMAKE_EXE%" -S . -B "%DFN_VS_DIR%" -G "Visual Studio 17 2022" -A x64 ^
    -DCMAKE_TOOLCHAIN_FILE="%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake" ^
    -DVCPKG_TARGET_TRIPLET=x64-windows ^
    -DVCPKG_OVERLAY_PORTS="%CD%\cmake\vcpkg-overlay" ^
    -DDFN_USE_KLU=ON
  goto :check
)

echo warning: SUNDIALS not found; generating a core-only solution.
echo          dfn-sim will compile but needs DFN_USE_KLU=ON to run.
"%CMAKE_EXE%" -S . -B "%DFN_VS_DIR%" -G "Visual Studio 17 2022" -A x64
goto :check

:build
"%CMAKE_EXE%" -S . -B "%DFN_VS_DIR%"

:check
if errorlevel 1 exit /b 1
if not exist "%DFN_VS_DIR%\numerical_p2d.sln" (
  echo error: solution was not generated
  exit /b 1
)

if defined DFN_VS_NO_OPEN (
  echo Solution ready: %CD%\%DFN_VS_DIR%\numerical_p2d.sln
  exit /b 0
)

echo Opening %CD%\%DFN_VS_DIR%\numerical_p2d.sln
start "" "%CD%\%DFN_VS_DIR%\numerical_p2d.sln"
