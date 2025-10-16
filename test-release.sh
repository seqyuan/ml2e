#!/bin/bash

# 测试Release配置的脚本
echo "🔍 检查Release配置..."

# 检查workflow文件是否存在
if [ -f ".github/workflows/release-static.yml" ]; then
    echo "✅ release-static.yml 存在"
else
    echo "❌ release-static.yml 不存在"
    exit 1
fi

# 检查是否使用了RELEASE_TOKEN
if grep -q "RELEASE_TOKEN" .github/workflows/release-static.yml; then
    echo "✅ 已配置使用 RELEASE_TOKEN"
else
    echo "❌ 未配置 RELEASE_TOKEN"
fi

# 检查文件路径是否正确
if grep -q "ml2e-linux-amd64/ml2e-linux-amd64" .github/workflows/release-static.yml; then
    echo "✅ 文件路径配置正确"
else
    echo "❌ 文件路径配置有问题"
fi

echo ""
echo "📋 下一步操作："
echo "1. 在GitHub仓库设置中添加 RELEASE_TOKEN secret"
echo "2. 推送新的tag来测试Release"
echo "3. 查看Actions日志确认是否成功"
