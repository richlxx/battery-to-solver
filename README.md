# 从一块电池到求解器

**写给零数值计算经验的你 · 以 P2D/DFN 电池模型为例**

这是一本"边跑代码边学数值计算"的开源书，也是一套能运行的 C++ 求解器。全书从一个四则运算就能算的 OCV 曲线出发，逐步走到完整的 DFN（P2D）模型、手写 BDF 与工业级 IDA + KLU 两套后端，最后用时间–电压曲线反推动力学与扩散参数。

不需要微积分证明，也不要求电化学基础；会 `if`、`for` 和函数调用即可。

## 你会得到什么

- 一本 12 章、每章带检查点、答案与可运行实验的教程。
- 一个可复现的 DFN 求解器：有限体积离散、解析 Jacobian、BDF2-DAE、事件驱动的 CC/CV、SPM 降阶与参数辨识。
- 一套验证方法：守恒检验、观测收敛阶、与独立实现的对比，以及参数可辨识性分析。

## 快速开始

Windows（需要 Visual Studio 的 C++ 工具链）：

```bat
scripts\build_msvc.bat
build\Release\ocv_curve.exe --n 101
```

Linux / macOS：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/ocv_curve --n 101
```

第一条曲线就能画出来：输出 `soc,sto_neg,sto_pos,ocv_v` 四列 CSV。

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
- [前言与阅读路线](book/README.md)
- 写作计划与里程碑：[PLAN.md](PLAN.md)

## 可选组件

- **IDA + KLU（第 8 章起）**：需要 vcpkg 与 SUNDIALS 的 KLU overlay，构建命令见 [附录 B](book/appendix/b-setup.md)，或直接运行 `scripts\build_msvc_klu.bat`。
- **CUDA 批量正演（第 11 章）**：需要 NVIDIA GPU 与 CUDA Toolkit，设置 `CUDA_PATH` 后运行 `scripts\build_cuda.bat`；没有 GPU 也能完成其余各章。

## 测试

```bat
ctest --test-dir build -C Release --output-on-failure
```

带 KLU 的构建把 `build` 换成 `build-klu`。书中的数值结论都有对应测试，可以亲手复现。

## 状态

本书与代码仍在编写中，欢迎通过 Issue 指出错误或提出改进建议。
