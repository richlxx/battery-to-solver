#!/usr/bin/env bash
# 编译 dfnFoam。优先使用系统的 openfoam-dev（$WM_PROJECT_DIR），
# 否则回退到仓库外的 OpenFOAM 源码头文件目录（通过 FOAM_SRC 指定）。
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

find_foam_src() {
    if [ -n "${FOAM_SRC:-}" ] && [ -d "$FOAM_SRC/src/OpenFOAM" ]; then
        echo "$FOAM_SRC"
        return
    fi
    if [ -n "${WM_PROJECT_DIR:-}" ] && [ -d "$WM_PROJECT_DIR/src/OpenFOAM" ]; then
        echo "$WM_PROJECT_DIR"
        return
    fi
    for candidate in \
        "$HOME"/openfoam-src/openfoam-OpenFOAM-* \
        "$HOME"/OpenFOAM/OpenFOAM-*/src/.. \
        /usr/lib/openfoam/openfoam*/src/..; do
        if [ -d "$candidate/src/OpenFOAM" ]; then
            echo "$candidate"
            return
        fi
    done
    return 1
}

FOAM_SRC_DIR="$(find_foam_src || true)"
if [ -z "$FOAM_SRC_DIR" ]; then
    echo "error: OpenFOAM headers not found." >&2
    echo "  install openfoam-dev, or set FOAM_SRC to a source tree" >&2
    echo "  (for this machine: FOAM_SRC=\$HOME/openfoam-src/openfoam-OpenFOAM-v1912.200626)" >&2
    exit 1
fi

# 用发布包的 wmakeLnInclude 生成头文件软链接目录；只做一次。
for dir in OpenFOAM finiteVolume Pstream meshTools OSspecific/POSIX; do
    if [ ! -d "$FOAM_SRC_DIR/src/$dir/lnInclude" ]; then
        "$FOAM_SRC_DIR/wmake/wmakeLnInclude" "$FOAM_SRC_DIR/src/$dir" >/dev/null
    fi
done

CXX_BIN="${CXX:-g++}"
if ! command -v "$CXX_BIN" >/dev/null 2>&1; then
    if [ -x "$HOME/toolchain/bin/g++" ]; then
        CXX_BIN="$HOME/toolchain/bin/g++"
    else
        echo "error: no C++ compiler found (set CXX)" >&2
        exit 1
    fi
fi

# Debian 的 OpenFOAM 只装运行库，MPI/zlib 可能没有 -dev 软链接，
# 这时直接链接带版本号的 .so。
LIB_PATH="/usr/lib:/usr/lib/x86_64-linux-gnu"
MPI_Z_LIBS=""
if [ -e /usr/lib/x86_64-linux-gnu/libmpi.so ] &&
   [ -e /usr/lib/x86_64-linux-gnu/libz.so ]; then
    MPI_Z_LIBS="-lmpi -lz"
else
    MPI_SO="$(ls /usr/lib/x86_64-linux-gnu/libmpi.so.* 2>/dev/null | head -1 || true)"
    Z_SO="$(ls /usr/lib/x86_64-linux-gnu/libz.so.* 2>/dev/null | head -1 || true)"
    if [ -z "$MPI_SO" ] || [ -z "$Z_SO" ]; then
        echo "error: MPI/zlib shared libraries not found" >&2
        exit 1
    fi
    MPI_Z_LIBS="$MPI_SO $Z_SO"
fi

INCLUDES=(
    -I"$FOAM_SRC_DIR/src/OpenFOAM/lnInclude"
    -I"$FOAM_SRC_DIR/src/OpenFOAM/include"
    -I"$FOAM_SRC_DIR/src/finiteVolume/lnInclude"
    -I"$FOAM_SRC_DIR/src/meshTools/lnInclude"
    -I"$FOAM_SRC_DIR/src/OSspecific/POSIX/lnInclude"
    -I"$FOAM_SRC_DIR/src/Pstream/lnInclude"
)

# shellcheck disable=SC2086
"$CXX_BIN" -std=c++14 -O2 -w \
    -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -DNoRepository \
    "${INCLUDES[@]}" "$here/dfnFoam.C" \
    -L/usr/lib -L/usr/lib/x86_64-linux-gnu \
    -lfiniteVolume -lmeshTools -lOpenFOAM -lPstream $MPI_Z_LIBS \
    -Wl,-rpath-link,"/usr/lib:/usr/lib/x86_64-linux-gnu" \
    -Wl,-rpath,/usr/lib/x86_64-linux-gnu \
    -ldl -lpthread -o "$here/dfnFoam"

echo "built $here/dfnFoam"
