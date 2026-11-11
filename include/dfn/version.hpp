#pragma once

#include <string_view>

namespace dfn {

// 库版本信息（语义化版本号），随发布与调试输出一起提供。
// 供命令行 --version 与各类报告头部使用。
inline constexpr std::string_view kVersion = "0.1.0";

}  // namespace dfn
