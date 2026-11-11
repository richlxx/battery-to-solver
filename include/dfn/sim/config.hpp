#pragma once

#include <string>

#include "dfn/params.hpp"

namespace dfn {

// 参数文件加载：把 configs/demo.yaml 的文本参数解析为 Params，
// 供命令行程序与实验共用同一套配置入口。

// 加载扁平 YAML 子集：顶层小节（cell、negative、positive、separator、
// electrolyte）内为 "key: number" 行。未知键会被拒绝，避免拼写错误静默
// 回退到默认值。
Params load_params_yaml(const std::string& path);

}  // namespace dfn
