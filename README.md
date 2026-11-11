# 从一块电池到求解器

**写给零数值计算经验的你 · 以 P2D/DFN 电池模型为例**

这是一本"边跑代码边学数值计算"的开源书，也是一套能运行的 C++ 求解器。全书从一个四则运算就能算的 OCV 曲线出发，逐步走到完整的 DFN（P2D）模型、手写 BDF 与工业级 IDA + KLU 两套后端，最后用时间–电压曲线反推动力学与扩散参数。

不需要微积分证明，也不要求电化学基础；会 `if`、`for` 和函数调用即可。

## 你会得到什么

- 一本 12 章、每章带检查点、答案与可运行实验的教程。
- 一个可复现的 DFN 求解器：有限体积离散、解析 Jacobian、BDF2-DAE、事件驱动的 CC/CV、SPM/SPMe/DFN 三级模型与参数辨识。
- 一套验证方法：守恒检验、观测收敛阶、与独立实现的对比，以及参数可辨识性分析。

## 怎么读

| 你的情况 | 建议 |
| --- | --- |
| 完全新手 | 从第 0 章顺序读，每章先运行代码再看公式 |
| 会编程但没学过数值方法 | 第 1 章快速浏览，重点第 2 章 |
| 有电池/电化学背景 | 第 3 章快读，重点第 4–9 章 |
| 只想先跑起来 | 读第 0.7 节，再跳到第 9 章 |

每章的代码都可以这样用：**先运行 → 改参数看变化 → 读关键函数 → 自己写一遍**。

## 快速开始

Windows（需要 Visual Studio 的 C++ 工具链）：

```bat
scripts\build_msvc.bat
build\Release\ocv_curve.exe --n 101
```

想在 Visual Studio 里浏览代码、断点调试，运行一次：

```bat
scripts\open_vs.bat
```

它会生成并打开 `build-vs\numerical_p2d.sln`。若已安装 SUNDIALS/KLU（`VCPKG_ROOT` 或本地 `vcpkg_installed`），脚本会自动启用 `DFN_USE_KLU`；否则生成不带 IDA/KLU 后端的解决方案，`dfn-sim` 仍可编译，但运行时会提示需要 KLU。

Linux / macOS：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/ocv_curve --n 101
```

也可以用一键脚本 `scripts/build_unix.sh`（加 `--test` 顺便跑测试）。

第一条曲线就能画出来：输出 `soc,sto_neg,sto_pos,ocv_v` 四列 CSV。

## 在 Visual Studio 中浏览与调试

解决方案按用途分成四个解决方案文件夹，避免几十个工程平铺在一起：

| 分组 | 内容 |
| --- | --- |
| Applications | `dfn-sim` 命令行 |
| Libraries | `dfn_core`、`dfn_ida`（启用 KLU 时） |
| Labs | 每章的单文件实验 |
| Tests | 单元测试与 CTest 辅助工程 |

`dfn_core` 内部再按 `Sources/models`、`Sources/solvers`、`Headers/models` 等过滤器组织，与仓库目录一一对应。头文件已加入工程，可以直接在解决方案资源管理器里跳转。

默认启动项目是 `dfn-sim`，按 **F5** 会用 Debug 配置运行 SPM 恒流放电示例，工作目录已设为仓库根目录，因此能直接读取 `configs/` 下的参数文件。改成 DFN/CCCV 等场景时，在“项目属性 → 调试 → 命令参数”里修改即可。要调试某个 lab，右键该实验选择“设为启动项目”后按 F5。

## 目录结构

```text
book/       全书正文（GitBook 源）与配图
labs/       每章的单文件实验
include/    求解器库头文件
src/        库实现与 dfn-sim 命令行
tests/      测试与验证代码
configs/    电池参数 YAML
scripts/    构建、画图与对比脚本
```

## 阅读入口

- [全书目录](book/SUMMARY.md)
- 写作计划与里程碑：[PLAN.md](PLAN.md)

## 可选组件

- **IDA + KLU（第 8 章起）**：Windows 需要 vcpkg 与 SUNDIALS 的 KLU overlay，运行 `scripts\build_msvc_klu.bat`；macOS 直接 `brew install cmake sundials` 后运行 `scripts/build_unix.sh --klu`。详见 [附录 B](book/appendix/b-setup.md)。
- **CUDA 批量正演（第 11 章）**：需要 NVIDIA GPU 与 CUDA Toolkit，设置 `CUDA_PATH` 后运行 `scripts\build_cuda.bat`；macOS（含 Apple Silicon）无 NVIDIA CUDA 支持，无法运行本章实验，但没有 GPU 也能完成其余各章。

## 测试

```bat
ctest --test-dir build -C Release --output-on-failure
```

带 KLU 的构建把 `build` 换成 `build-klu`。书中的数值结论都有对应测试，可以亲手复现。

## 状态

本书与代码仍在编写中，欢迎通过 Issue 指出错误或提出改进建议。
