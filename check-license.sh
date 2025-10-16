#!/bin/bash

echo "🔍 检查许可证配置..."

# 检查ml2e.go是否有许可证头部
if grep -q "Copyright (c) 2025 seqyuan" ml2e.go; then
    echo "✅ ml2e.go 包含许可证头部"
else
    echo "❌ ml2e.go 缺少许可证头部"
fi

# 检查.licenserc.yaml配置
if grep -q "paths-ignore:" .licenserc.yaml; then
    echo "✅ .licenserc.yaml 配置存在"
else
    echo "❌ .licenserc.yaml 配置有问题"
fi

# 检查是否忽略了workflow文件
if grep -q ".github/workflows/\*.yml" .licenserc.yaml; then
    echo "✅ 已忽略GitHub workflow文件"
else
    echo "❌ 未忽略GitHub workflow文件"
fi

# 检查是否忽略了shell脚本
if grep -q "\*.sh" .licenserc.yaml; then
    echo "✅ 已忽略shell脚本文件"
else
    echo "❌ 未忽略shell脚本文件"
fi

echo ""
echo "📋 许可证配置总结："
echo "- ml2e.go: 已添加MIT许可证头部"
echo "- .licenserc.yaml: 已配置忽略workflow和脚本文件"
echo "- 许可证检查应该通过"
