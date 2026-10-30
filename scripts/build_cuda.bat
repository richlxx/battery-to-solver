@echo off
setlocal

cd /d "%~dp0.."

set "DFN_BUILD_DIR=build-cuda"
call "%~dp0build_msvc.bat" -DDFN_USE_CUDA=ON %*
