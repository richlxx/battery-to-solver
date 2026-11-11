// pybind11 绑定入口：把求解器的少量接口暴露成 Python 扩展模块 dfn_py。
// 目前只导出编译期版本号，供 Python 侧确认加载的扩展与源码版本一致。

#include <string>

#include <pybind11/pybind11.h>

#include "dfn/version.hpp"

// 定义 Python 扩展模块 dfn_py。m 为该模块对象，后续的类与函数绑定
// 都应写在这个函数体内。
PYBIND11_MODULE(dfn_py, m) {
  // 把 version.hpp 中的 kVersion 暴露为模块属性 __version__（Python str）。
  m.attr("__version__") = std::string(dfn::kVersion);
}
