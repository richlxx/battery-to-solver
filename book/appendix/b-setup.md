# 附录 B 环境、构建与 C++ 最小入门

## B.1 安装清单

| 组件 | 用途 | 获取方式 |
| --- | --- | --- |
| CMake ≥ 3.24 | 构建 | 官网安装包 |
| C++20 编译器 | 编译 | MSVC 2022 / GCC 12+ / Clang 15+ |
| vcpkg | 依赖管理 | 官方 bootstrap 脚本 |
| Eigen3 | 线性代数 | vcpkg |
| SUNDIALS（KLU 特性） | IDA + KLU | vcpkg `sundials[klu]` |
| pybind11、yaml-cpp、Catch2 | 绑定/配置/测试 | vcpkg |
| uv | Python 环境 | 独立安装 |

## B.2 第一个程序：不装任何依赖

第 0–2 章的 lab 只依赖标准库，可以直接编译：

```bash
g++ -std=c++20 labs/ocv_curve.cpp -o ocv
./ocv > ocv.csv
python scripts/plot_ocv.py ocv.csv
```

Windows 上如果只有 MSVC，在“Developer Command Prompt”里用：

```bat
cl /std:c++20 /EHsc labs\ocv_curve.cpp
ocv_curve.exe > ocv.csv
```

## B.3 完整项目构建

Windows 一键构建（推荐，使用 Visual Studio 自带的 CMake 4.x）：

```bat
scripts\build_msvc.bat
```

脚本会通过 `vswhere` 找到 VS、使用 VS 自带的 CMake 并选择对应的 Visual Studio 生成器。**CMake ≤ 3.31 不认识 VS 18，且在中文路径下可能直接崩溃**，所以不要用系统旧版 CMake 配置本仓库。

### 带 IDA + KLU 的构建（第 8 章）

SUNDIALS 的官方 vcpkg 端口不包含 KLU 特性，本仓库提供了一份 overlay 端口，把 `klu` 映射到 `SUNDIALS_ENABLE_KLU`：

```bat
vcpkg install --classic suitesparse-klu:x64-windows
vcpkg install --classic --overlay-ports=cmake/vcpkg-overlay "sundials[klu]:x64-windows"
scripts\build_msvc_klu.bat
ctest --test-dir build-klu -C Release --output-on-failure
```

- `suitesparse-klu` 会连带编译 OpenBLAS 与 AMD/COLAMD/BTF，第一次安装约 4–5 分钟。
- SUNDIALS 默认按动态库构建；CMake 会在构建后把用到的 DLL 复制到可执行文件旁边，所以运行测试不需要改 `PATH`。
- `build_msvc_klu.bat` 会自动打开 `-DDFN_USE_KLU=ON` 并把 `CMAKE_PREFIX_PATH` 指向 vcpkg 的 `installed/x64-windows`。
- 若配置时报 `SUNDIALS::sunlinsolklu is missing`，说明 vcpkg 装的是不带 KLU 的官方包，先卸载再按上面的 overlay 命令重装。

手动构建（Linux/macOS，或英文路径下的 Windows）：

```bash
cmake -S . -B build
cmake --build build --config Release -j
ctest --test-dir build -C Release --output-on-failure
```

带 KLU 时加上：

```bash
cmake -S . -B build-klu -DDFN_USE_KLU=ON \
  -DCMAKE_PREFIX_PATH="$VCPKG_ROOT/installed/x64-windows"
```

注意：不装 vcpkg 也能**构建**全部章节的代码与测试；第 8 章之后的仿真**运行**需要 SUNDIALS（见下）。参数文件用的是自带的 YAML 子集解析器，不依赖第三方库。

## B.4 Python 侧

第三层验证需要 PyBaMM：

```bash
uv venv .venv
uv pip install -i https://pypi.tuna.tsinghua.edu.cn/simple pybamm matplotlib
```

先让 C++ 程序生成曲线，再对比：

```bash
build-klu/Release/dfn-sim --mode cc --crate 1 --dt 60 --out results/cc.csv
.venv/Scripts/python.exe scripts/compare_pybamm.py results/cc.csv --save results/cc.png
```

PyBaMM 首次导入较慢，属正常；参数映射写在 `scripts/compare_pybamm.py` 里，改参数时两边要一起改。

## B.5 C++ 最小入门

只列书中会用到的语法，够读代码即可。

```cpp
#include <iostream>
#include <vector>

// 函数：输入半径，返回球体积
double sphere_volume(double r) {
  return 4.0 / 3.0 * 3.141592653589793 * r * r * r;
}

int main() {
  std::vector<double> radii = {1e-6, 2e-6, 5e-6};
  for (double r : radii) {
    std::cout << r << " -> " << sphere_volume(r) << "\n";
  }
}
```

| 语法 | 含义 | 书中出现处 |
| --- | --- | --- |
| `std::vector<double>` | 可变长数组 | 存浓度、网格 |
| `for (double x : v)` | 遍历 | 装配循环 |
| `double f(double)` | 函数 | 算子、残差 |
| `struct Params { ... };` | 打包一组数据 | 参数集 |
| `const Params& p` | 只读引用，避免拷贝 | 函数参数 |
| `std::ofstream` | 写文件 | 输出 CSV |
| `throw std::runtime_error("...")` | 报错 | 参数校验 |

读编译错误的顺序：先看第一条错误 → 找到行号 → 检查括号与分号 → 再看类型是否匹配。后面的错误常常是第一条的连锁反应。

## B.6 无 KLU 的降级构建

```bash
cmake -S . -B build-noklu -DDFN_USE_KLU=OFF
```

这是默认配置。第 8 章的 `dfn_ida`、`klu_probe`、`ida_compare` 与 `dfn_ida_tests` 不会生成，其余章节的代码与测试不受影响。`dfn-sim` 仍可构建，但运行仿真时会提示"built without SUNDIALS"；参数文件加载和 `--help` 不受影响。第 8、9 章的正文仍可阅读，只是不能运行 IDA 相关实验。

## B.7 常见问题

- SUNDIALS 版本与 `SUNContext` API 不匹配：本项目按 7.8.0 编写；升级前先查发布说明，注意 `IDAInit`、`IDASetId` 等签名变化。
- pybind11 与 PyBaMM 的 Python 版本冲突：两者使用独立 venv。
- SuiteSparse 链接失败：确认 vcpkg 三元组与编译器架构一致（x64-windows）。
- 中文路径导致构建失败：把仓库放在纯英文路径下，或使用 VS 自带的 CMake 4.x。
- 运行时报 `找不到 sundials_ida.dll`：检查 DLL 是否被复制到 `build-klu\Release`，或把 `%VCPKG_ROOT%\installed\x64-windows\bin` 加入 `PATH`。
