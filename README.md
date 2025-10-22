# ml2e - FastQ read ID modifier

一个基于Go的高性能FastQ文件readID修改工具，将readID中的`@ML15`前缀替换为`@E251`前缀。

## 功能特性

- **高性能批量处理**：支持多线程并行处理，充分利用CPU资源
- **压缩文件支持**：原生支持gzip压缩文件格式
- **动态缓存机制**：根据worker数量自动调整缓存大小，优化内存使用
- **pigz集成**：支持使用pigz进行快速压缩和解压缩

## 使用方法

### 基本语法

```bash
ml2e [选项]
```

### 必需参数

- `-fq1`：输入FastQ文件1（支持`.gz`格式）
- `-fq2`：输入FastQ文件2（支持`.gz`格式）
- `-outdir`：输出目录

### 可选参数

- `-workers`：worker线程数量（默认：4，建议6-8）
- `-pigz`：pigz可执行文件路径，用于压缩输出文件
- `-pigz-decompress`：使用pigz进行解压缩（对gzip文件更快）
- `-version`：显示版本信息

### 使用示例

#### 基本使用
```bash
ml2e -fq1 input1.fastq.gz -fq2 input2.fastq.gz -outdir output
```

#### 高并发处理
```bash
ml2e -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output -workers 8
```

#### 使用pigz压缩输出
```bash
ml2e -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output -pigz /usr/bin/pigz
```

#### 使用pigz解压缩
```bash
ml2e -fq1 input1.fastq.gz -fq2 input2.fastq.gz -outdir output -pigz /usr/bin/pigz -pigz-decompress
```

## 输出文件

程序会在指定的输出目录中生成以下文件：

1. **修改后的FastQ文件**：
   - `{原文件名1}`：修改后的R1文件
   - `{原文件名2}`：修改后的R2文件
   - 如果指定了`-pigz`参数，文件将被压缩为`.gz`格式

2. **统计报告**：
   - `{原文件名2}.stats.txt`：包含处理统计信息

### 统计报告内容

统计报告包含以下信息：
- 总reads处理数量
- 总输出reads数量

## 工作原理

程序会将所有readID中开头的`@ML15`替换为`@E251`，其他内容保持不变。例如：
- `@ML1512345:1:1101:1000:2211/1` → `@E25112345:1:1101:1000:2211/1`
- `@AB12345:1:1101:1000:2211/1` → `@AB12345:1:1101:1000:2211/1`（保持不变）

所有reads都会被保留并输出，只是readID会被修改。


## 安装

### 使用 go install 安装

```bash
# 安装最新版本
go install github.com/seqyuan/ml2e@latest

# 安装特定版本
go install github.com/seqyuan/ml2e@v0.1.7
```

安装后，程序会位于 `$GOPATH/bin/ml2e`，确保该目录在您的PATH环境变量中。

### 从源码构建

```bash
git clone https://github.com/seqyuan/ml2e.git
cd ml2e
go build -o ml2e ml2e.go
```

## Release

```shell
git add -A
git commit -m "release"
git tag v0.1.7
git push origin main
git push origin v0.1.7
```