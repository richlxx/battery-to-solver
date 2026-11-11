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

rem The generator name must match the installed IDE, otherwise MSBuild looks
rem for the old toolset (e.g. v143 under VS 18) and fails with MSB8020.
for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationVersion`) do set "VSVER=%%i"
if not defined VSVER (
  echo error: cannot determine Visual Studio version
  exit /b 1
)
for /f "delims=. tokens=1" %%i in ("%VSVER%") do set "VSMAJOR=%%i"

set "VSGEN="
if "%VSMAJOR%"=="18" set "VSGEN=Visual Studio 18 2026"
if "%VSMAJOR%"=="17" set "VSGEN=Visual Studio 17 2022"
if "%VSMAJOR%"=="16" set "VSGEN=Visual Studio 16 2019"
if not defined VSGEN (
  echo error: unsupported Visual Studio version %VSVER%
  exit /b 1
)

set "FRESH="
set "CMAKE_EXE=%VSPATH%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if not exist "%CMAKE_EXE%" (
  echo error: CMake bundled with Visual Studio not found
  exit /b 1
)

if not defined DFN_VS_DIR set "DFN_VS_DIR=build-vs"

set "CACHED_GEN="
if exist "%DFN_VS_DIR%\CMakeCache.txt" for /f "tokens=2 delims==" %%i in ('findstr /b /c:"CMAKE_GENERATOR:INTERNAL=" "%DFN_VS_DIR%\CMakeCache.txt"') do set "CACHED_GEN=%%i"
if defined CACHED_GEN if not "%CACHED_GEN%"=="%VSGEN%" (
  echo Build cache uses "%CACHED_GEN%", rebuilding with "%VSGEN%".
  set "FRESH=--fresh"
)

if exist "%DFN_VS_DIR%\CMakeCache.txt" if not defined FRESH goto :build

rem First configure: enable the SUNDIALS/KLU backend when it is already
rem installed, so that dfn-sim (the startup project) can actually run.
set "SUNDIALS_PREFIX="
if defined VCPKG_ROOT if exist "%VCPKG_ROOT%\installed\x64-windows\share\sundials\SUNDIALSConfig.cmake" set "SUNDIALS_PREFIX=%VCPKG_ROOT%\installed\x64-windows"
if not defined SUNDIALS_PREFIX if exist "%CD%\vcpkg_installed\x64-windows\share\sundials\SUNDIALSConfig.cmake" set "SUNDIALS_PREFIX=%CD%\vcpkg_installed\x64-windows"
if not defined SUNDIALS_PREFIX if exist "%CD%\build-klu\vcpkg_installed\x64-windows\share\sundials\SUNDIALSConfig.cmake" set "SUNDIALS_PREFIX=%CD%\build-klu\vcpkg_installed\x64-windows"

if defined SUNDIALS_PREFIX (
  echo Configuring with SUNDIALS/KLU: %SUNDIALS_PREFIX%
  "%CMAKE_EXE%" -S . -B "%DFN_VS_DIR%" %FRESH% -G "%VSGEN%" -A x64 ^
    -DDFN_USE_KLU=ON -DCMAKE_PREFIX_PATH="%SUNDIALS_PREFIX%"
  goto :check
)

if defined VCPKG_ROOT (
  echo SUNDIALS is not installed yet; letting vcpkg install the manifest.
  "%CMAKE_EXE%" -S . -B "%DFN_VS_DIR%" %FRESH% -G "%VSGEN%" -A x64 ^
    -DCMAKE_TOOLCHAIN_FILE="%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake" ^
    -DVCPKG_TARGET_TRIPLET=x64-windows ^
    -DVCPKG_OVERLAY_PORTS="%CD%\cmake\vcpkg-overlay" ^
    -DDFN_USE_KLU=ON
  goto :check
)

echo warning: SUNDIALS not found; generating a core-only solution.
echo          dfn-sim will compile but needs DFN_USE_KLU=ON to run.
"%CMAKE_EXE%" -S . -B "%DFN_VS_DIR%" %FRESH% -G "%VSGEN%" -A x64
goto :check

:build
"%CMAKE_EXE%" -S . -B "%DFN_VS_DIR%" %FRESH% -G "%VSGEN%" -A x64

:check
if errorlevel 1 exit /b 1
rem VS 18 writes the new .slnx format; older IDEs write .sln.
set "SLN=%DFN_VS_DIR%\numerical_p2d.sln"
if not exist "%SLN%" set "SLN=%DFN_VS_DIR%\numerical_p2d.slnx"
if not exist "%SLN%" (
  echo error: solution was not generated
  exit /b 1
)

if defined DFN_VS_NO_OPEN (
  echo Solution ready: %CD%\%SLN%
  exit /b 0
)

echo Opening %CD%\%SLN%
if exist "%VSPATH%\Common7\IDE\devenv.exe" (
  start "" "%VSPATH%\Common7\IDE\devenv.exe" "%CD%\%SLN%"
) else (
  start "" "%CD%\%SLN%"
)
