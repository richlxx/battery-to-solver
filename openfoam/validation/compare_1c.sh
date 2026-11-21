#!/usr/bin/env bash
# 运行 dfnFoam 的 1C 算例，并与仓库根目录的 dfn-sim 参考解逐时刻比较。
# 参考解需要事先构建（例如 build-ref/dfn-sim），也可以用 DFN_SIM 指定。
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "${here}/../.." && pwd)"
case_dir="${root}/openfoam/cases/demo-1c"
dfn_sim="${DFN_SIM:-${root}/build-ref/dfn-sim}"
ref_csv="${TMPDIR:-/tmp}/dfn_sim_1c.csv"
of_csv="${case_dir}/postProcessing/dfnFoam/terminalVoltage.csv"

if [[ ! -x "${dfn_sim}" ]]; then
  echo "error: reference solver not found: ${dfn_sim}" >&2
  echo "       build it first or set DFN_SIM=/path/to/dfn-sim" >&2
  exit 1
fi

"${case_dir}/run.sh" >"${TMPDIR:-/tmp}/dfnFoam_1c.log" 2>&1
"${dfn_sim}" --model dfn --mode cc --crate 1 --t-max 3600 --dt 600 \
  --nx 8 4 8 --nr 6 6 --out "${ref_csv}" --quiet 2>/dev/null

printf '%6s %11s %11s %9s %10s %10s\n' \
  "t(s)" "dfnFoam V" "dfn-sim V" "dV(mV)" "theta_neg" "theta_pos"
awk -F, '
  NR == FNR {
    if (FNR > 1 && $1 > 0 && $1 % 600 == 0) {
      rv[$1] = $2
    }
    next
  }
  FNR > 1 && ($1 in rv) {
    printf "%6.0f %11.5f %11.5f %+9.2f %10.5f %10.5f\n",
           $1, $2, rv[$1], 1000.0 * ($2 - rv[$1]), $4, $5
  }
' "${ref_csv}" "${of_csv}"

echo "参考解记录: ${ref_csv}"
echo "dfnFoam 全量结果: ${of_csv}"
