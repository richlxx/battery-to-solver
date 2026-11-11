// 测试文件：构建冒烟测试（对应全局 version 模块）。
// 验证测试工程能够链接 dfn 库，且导出的版本字符串符合预期。

#include <iostream>
#include <string>

#include "dfn/version.hpp"

int main() {
  // 版本号是 ABI/配置兼容性的第一道检查；与发布约定 0.1.0 不一致即失败。
  const std::string version(dfn::kVersion);
  if (version != "0.1.0") {
    std::cerr << "unexpected version: " << version << "\n";
    return 1;
  }
  std::cout << "dfn smoke ok (" << version << ")\n";
  return 0;
}
