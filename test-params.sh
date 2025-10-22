#!/bin/bash

echo "🧪 测试新的参数结构..."

# 测试参数验证
echo "📋 测试参数验证："
./ml2e 2>&1 | head -1

echo ""
echo "📋 测试帮助信息："
./ml2e 2>&1 | grep -E "(-fq1|-fq2|-r1|-r2|-from|-to)"

echo ""
echo "✅ 新参数结构测试完成！"
echo ""
echo "📋 使用方法示例："
echo "ml2e -fq1 input1.fastq.gz -fq2 input2.fastq.gz -r1 output1.fastq -r2 output2.fastq"
echo "ml2e -fq1 input1.fastq -fq2 input2.fastq -r1 output1.fastq -r2 output2.fastq -from ML15 -to E251"
