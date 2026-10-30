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

rem Use the CMake bundled with Visual Studio: newer versions know the current
rem MSVC toolset and handle non-ASCII source paths.
set "CMAKE_EXE=%VSPATH%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if not exist "%CMAKE_EXE%" (
  echo error: CMake bundled with Visual Studio not found
  exit /b 1
)

if not defined DFN_BUILD_DIR set "DFN_BUILD_DIR=build"

"%CMAKE_EXE%" -S . -B "%DFN_BUILD_DIR%" %*
if errorlevel 1 exit /b 1

"%CMAKE_EXE%" --build "%DFN_BUILD_DIR%" --config Release
