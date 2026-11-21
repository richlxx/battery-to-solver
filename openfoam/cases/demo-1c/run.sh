#!/usr/bin/env bash
# 构建（若需要）并运行 1C 恒流放电算例。
# 需要 blockMesh 在 PATH 中；Debian 的 OpenFOAM 安装还需要 WM_PROJECT_DIR。
set -euo pipefail

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
solver_dir="$(cd "${case_dir}/../../dfnFoam" && pwd)"
export WM_PROJECT_DIR="${WM_PROJECT_DIR:-/usr/share/openfoam}"

if [[ ! -x "${solver_dir}/dfnFoam" ]]; then
  "${solver_dir}/build.sh"
fi

cd "${case_dir}"
# 删除上一次运行留下的时间目录（保留 0），避免 startFrom startTime 读到旧场。
find . -maxdepth 1 -type d -regex '\./[1-9][0-9]*' -exec rm -rf {} +
rm -rf postProcessing

blockMesh > log.blockMesh 2>&1
"${solver_dir}/dfnFoam" | tee log.dfnFoam

echo "结果: ${case_dir}/postProcessing/dfnFoam/terminalVoltage.csv"
