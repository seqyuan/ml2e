# PatternQC - FASTQ 过滤含有固定序列reads的工具

一个基于go的高性能 FASTQ 文件过滤工具，支持多线程处理和压缩文件格式。

## 安装

### 方法一：使用 go install（推荐）

```bash
# 安装最新版本
go install github.com/seqyuan/patternqc@latest

# 安装特定版本
go install github.com/seqyuan/patternqc@v0.8.0
```

安装后，`patternqc` 命令会被安装到 `$GOPATH/bin` 目录中，确保该目录在你的 `PATH` 环境变量中。

### 方法二：从源码编译

#### 依赖要求

- Go 1.22 或更高版本
- `github.com/seqyuan/annogene` 包

#### 编译步骤

```bash
# 克隆仓库
git clone https://github.com/seqyuan/patternqc.git
cd patternqc

# 编译
go build -o patternqc patternqc.go
```

### 验证安装

```bash
# 检查版本
patternqc -version

# 查看帮助
patternqc
```

## 使用方法

### 基本语法

```bash
./patternqc [选项]
```

### 必需参数

- `-fq1`：输入 FASTQ 文件 1（支持 `.gz` 格式）
- `-fq2`：输入 FASTQ 文件 2（支持 `.gz` 格式）
- `-outdir`：输出目录

### 可选参数

- `-pattern`：要搜索的 pattern（默认：`AGCAGTGGTATCAACGCAGAGTACA`）
- `-percent`：保留 pattern reads 的百分比（0-100，默认：5）
- `-workers`：worker线程数量（默认：4，建议6-8）
- `-pigz`：pigz 可执行文件路径，用于压缩输出文件

### 使用示例

#### 基本使用
```bash
patternqc -fq1 ./data/f1.fq.gz -fq2 ./data/f2.fq.gz -outdir ./result
```

#### 自定义参数
```bash
patternqc -fq1 ./data/f1.fq.gz -fq2 ./data/f2.fq.gz -outdir ./result -percent 10
```

#### 自定义 pattern
```bash
patternqc -fq1 ./data/f1.fq.gz -fq2 ./data/f2.fq.gz -outdir ./result -pattern "AGCAGTGGTATCAACGCAGAGTACA" -percent 5
```

#### 使用 pigz 压缩输出
```bash
patternqc -fq1 ./data/f1.fq.gz -fq2 ./data/f2.fq.gz -outdir ./result -pigz /usr/bin/pigz
```

#### 高并发处理
```bash
patternqc -fq1 ./data/f1.fq.gz -fq2 ./data/f2.fq.gz -outdir ./result -workers 8
```

#### 自定义pattern和高并发
```bash
patternqc -fq1 ./data/f1.fq.gz -fq2 ./data/f2.fq.gz -outdir ./result -pattern "AGCAGTGGTATCAACGCAGAGTACA" -percent 10 -workers 6
```

## 输出文件

程序会在指定的输出目录中生成以下文件：

1. **过滤后的 FASTQ 文件**：
   - `{原文件名1}`：过滤后的 R1 文件
   - `{原文件名2}`：过滤后的 R2 文件
   - 如果指定了 `-pigz` 参数，文件将被压缩为 `.gz` 格式

2. **统计报告**：
   - `{原文件名2}.stats.txt`：包含处理统计信息

### 统计报告内容

统计报告包含以下信息：
- 总 reads 数量
- 原始 pattern reads 数量
- 保留的 pattern reads 数量
- 最终比例


### 自动发布流程

1. **创建新版本标签**：
   ```bash
   git tag v0.8.0
   git push origin v0.8.0
   ```

2. **自动触发Release**：
   - 推送tag后自动触发GitHub Actions
   - 构建多平台二进制文件（Linux, Windows, macOS）
   - 自动创建GitHub Release并上传文件

## 技术支持

如有问题或建议，请联系开发团队。
