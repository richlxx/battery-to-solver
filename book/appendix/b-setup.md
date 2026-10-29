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

手动构建（Linux/macOS，或英文路径下的 Windows）：

```bash
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
cmake --build build --config Release -j
ctest --test-dir build -C Release --output-on-failure
```

注意：M0 阶段的所有目标（`dfn-sim`、`ocv_curve`、冒烟测试）**只依赖标准库**，不装 vcpkg 也能构建。Eigen、YAML、SUNDIALS 从 M2/M4 才开始使用。

## B.4 Python 侧

```bash
uv venv .venv
uv pip install numpy matplotlib pybamm
uv run python scripts/compare_pybamm.py --config configs/chen2020.yaml
```

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

第 5、8 章改用 banded 线性求解器，功能完整但性能下降；书中给出对比数据。

## B.7 常见问题

- SUNDIALS 版本与 `SUNContext` API 不匹配：固定 7.1.x，升级前先查发布说明。
- pybind11 与 PyBaMM 的 Python 版本冲突：两者使用独立 venv。
- SuiteSparse 链接失败：确认 vcpkg 三元组与编译器架构一致（x64-windows）。
- 中文路径导致构建失败：把仓库放在纯英文路径下。
