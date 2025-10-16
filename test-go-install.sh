#!/bin/bash

echo "🔍 测试 go install 功能..."

# 检查当前目录结构
echo "📁 当前目录结构："
ls -la

echo ""
echo "📋 检查go.mod配置："
cat go.mod

echo ""
echo "🔧 测试构建："
if go build -v .; then
    echo "✅ 本地构建成功"
else
    echo "❌ 本地构建失败"
    exit 1
fi

echo ""
echo "📦 检查二进制文件："
if [ -f "ml2e" ]; then
    echo "✅ 二进制文件存在"
    ./ml2e -version
else
    echo "❌ 二进制文件不存在"
fi

echo ""
echo "🚀 go install 命令示例："
echo "go install github.com/seqyuan/ml2e@latest"
echo "go install github.com/seqyuan/ml2e@v0.1.3"
echo ""
echo "📋 前提条件："
echo "1. 确保仓库是公开的"
echo "2. 确保已创建并推送了v0.1.3标签"
echo "3. 确保go.mod文件配置正确"
echo "4. 确保有main函数"
