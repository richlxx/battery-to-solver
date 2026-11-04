#!/usr/bin/env bash
# Unix (Linux / macOS) 一键构建，对应 Windows 的 build_msvc.bat。
#
# 用法：
#   scripts/build_unix.sh                # 基础构建到 build/
#   scripts/build_unix.sh --klu          # IDA + KLU 构建到 build-klu/
#   scripts/build_unix.sh --dir mybuild  # 自定义构建目录
#   scripts/build_unix.sh --test         # 构建后运行 ctest
# 其余参数原样传给 cmake 配置（如 -DDFN_USE_KLU=ON 会自动改用 build-klu/）。
#
# macOS 上 KLU 后端直接使用 Homebrew 的 sundials（自带 KLU 支持）：
#   brew install cmake sundials
set -euo pipefail

cd "$(dirname "$0")/.."

build_dir="build"
args=()
run_tests=0
while [ $# -gt 0 ]; do
  case "$1" in
    --klu)
      build_dir="build-klu"
      args+=(-DDFN_USE_KLU=ON)
      ;;
    --dir)
      build_dir="$2"
      shift
      ;;
    --test)
      run_tests=1
      ;;
    *)
      if [ "$1" = "-DDFN_USE_KLU=ON" ]; then
        build_dir="build-klu"
      fi
      args+=("$1")
      ;;
  esac
  shift
done

if [ "$(uname -s)" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
  # Homebrew 的 sundials 公式自带 KLU；显式指向前缀，避免 PATH 不含 brew 时找不到。
  args+=("-DCMAKE_PREFIX_PATH=$(brew --prefix)")
  if [ "$build_dir" = "build-klu" ] && [ -d "$(brew --prefix)/opt/libomp" ]; then
    # AppleClang 不带 OpenMP，指向 Homebrew 的 libomp（sundials 安装时已作为依赖装好）。
    args+=("-DOpenMP_ROOT=$(brew --prefix)/opt/libomp")
  fi
fi

cmake -S . -B "$build_dir" "${args[@]}"
cmake --build "$build_dir" -j

if [ "$run_tests" = 1 ]; then
  ctest --test-dir "$build_dir" --output-on-failure
fi

echo "构建完成：$build_dir"
