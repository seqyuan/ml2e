package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/seqyuan/annogene/io/fastq"
)

const (
	Version   = "0.1.2"
	BuildTime = "2025-10-16"
)

// 简化的结果结构 - 只需要配对信息
type ProcessResult struct {
	Seq1 fastq.Sequence
	Seq2 fastq.Sequence
}

// 修改readID中的@ML为@E
func modifyReadID(readID string) string {
	if strings.HasPrefix(readID, "@ML") {
		return "@E" + readID[3:] // 将@ML替换为@E
	}
	return readID
}

// 获取文件basename（保持原始扩展名，但移除.gz）
func getBasename(filePath string) string {
	filename := filepath.Base(filePath)

	// 移除.gz扩展名
	if strings.HasSuffix(filename, ".gz") {
		filename = strings.TrimSuffix(filename, ".gz")
	}

	return filename
}

// 创建reader，支持普通文件和gz文件
func createReader(filePath string) (io.Reader, *os.File, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to open file %s: %v", filePath, err)
	}

	var reader io.Reader = file

	// 检查是否为gz文件
	if strings.HasSuffix(filePath, ".gz") {
		// 使用gzip reader，增加缓冲区大小
		gz, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, nil, fmt.Errorf("failed to create gzip reader for %s: %v", filePath, err)
		}
		// 使用更大的缓冲区包装gzip reader
		reader = bufio.NewReaderSize(gz, 1024*1024) // 1MB缓冲区
	}

	return reader, file, nil
}

// 创建writer - 只输出文本文件，不压缩
func createWriter(filePath string) (io.Writer, *os.File, error) {
	// 创建输出目录
	dir := filepath.Dir(filePath)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return nil, nil, fmt.Errorf("failed to create output directory: %v", err)
	}

	file, err := os.Create(filePath)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to create output file %s: %v", filePath, err)
	}

	return file, file, nil
}

// 处理单个read pair的worker函数 - 修改readID
func processReadPair(seq1, seq2 fastq.Sequence) *ProcessResult {
	// 修改seq1的readID
	modifiedSeq1 := fastq.Sequence{
		ID1:     []byte(modifyReadID(string(seq1.ID1))),
		Letters: seq1.Letters,
		ID2:     seq1.ID2,
		Quality: seq1.Quality,
	}

	// 修改seq2的readID
	modifiedSeq2 := fastq.Sequence{
		ID1:     []byte(modifyReadID(string(seq2.ID1))),
		Letters: seq2.Letters,
		ID2:     seq2.ID2,
		Quality: seq2.Quality,
	}

	return &ProcessResult{
		Seq1: modifiedSeq1,
		Seq2: modifiedSeq2,
	}
}

// 批量处理reads的worker函数
func processReadPairBatch(readPairs []ReadPair) []*ProcessResult {
	results := make([]*ProcessResult, 0, len(readPairs))

	for _, readPair := range readPairs {
		result := processReadPair(readPair.Seq1, readPair.Seq2)
		results = append(results, result)
	}

	return results
}

// 使用pigz压缩文件
func compressWithPigz(file1, file2, pigzPath string, numWorkers int) error {
	// 检查pigz路径是否存在
	if _, err := os.Stat(pigzPath); os.IsNotExist(err) {
		return fmt.Errorf("pigz not found at path: %s", pigzPath)
	}

	// 使用goroutine同时压缩两个文件
	var wg sync.WaitGroup
	var err1, err2 error

	wg.Add(2)

	// 启动第一个文件的压缩
	go func() {
		defer wg.Done()
		err1 = compressFileWithPigz(file1, pigzPath, numWorkers)
		if err1 != nil {
			err1 = fmt.Errorf("failed to compress %s: %v", file1, err1)
		}
	}()

	// 启动第二个文件的压缩
	go func() {
		defer wg.Done()
		err2 = compressFileWithPigz(file2, pigzPath, numWorkers)
		if err2 != nil {
			err2 = fmt.Errorf("failed to compress %s: %v", file2, err2)
		}
	}()

	// 等待两个压缩任务完成
	wg.Wait()

	// 检查是否有错误发生
	if err1 != nil {
		return err1
	}
	if err2 != nil {
		return err2
	}

	return nil
}

// 使用pigz解压缩gzip文件
func decompressWithPigz(filePath, pigzPath string, numWorkers int) (string, error) {
	// 检查pigz路径是否存在
	if _, err := os.Stat(pigzPath); os.IsNotExist(err) {
		return filePath, fmt.Errorf("pigz not found at path: %s", pigzPath)
	}

	// 创建临时文件
	tempFile, err := os.CreateTemp("", "ml2e_*.fq")
	if err != nil {
		return filePath, fmt.Errorf("failed to create temp file: %v", err)
	}
	tempFilePath := tempFile.Name()
	tempFile.Close()

	// 构建pigz命令
	cmd := exec.Command(pigzPath, "-d", "-p", fmt.Sprintf("%d", numWorkers), "-o", tempFilePath, filePath)

	// 设置输出
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	// 执行解压缩
	if err := cmd.Run(); err != nil {
		os.Remove(tempFilePath) // 清理临时文件
		return filePath, fmt.Errorf("failed to decompress %s: %v", filePath, err)
	}

	fmt.Printf("Decompressed %s to %s using pigz\n", filePath, tempFilePath)
	return tempFilePath, nil
}

// 使用pigz压缩单个文件
func compressFileWithPigz(filePath, pigzPath string, numWorkers int) error {
	// 构建pigz命令
	cmd := exec.Command(pigzPath, "-p", fmt.Sprintf("%d", numWorkers), filePath)

	// 设置输出
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	// 执行命令
	return cmd.Run()
}

// 高性能缓冲写入器 - 优化锁粒度
type BufferedFastqWriter struct {
	mu      sync.Mutex
	buffer  []byte
	writer  fastq.Writer
	file    *os.File
	maxSize int
	written int64 // 使用int64支持原子操作
}

func NewBufferedFastqWriter(writer fastq.Writer, file *os.File, bufferSize int) *BufferedFastqWriter {
	return &BufferedFastqWriter{
		buffer:  make([]byte, 0, bufferSize),
		writer:  writer,
		file:    file,
		maxSize: bufferSize,
	}
}

func (bfw *BufferedFastqWriter) Write(seq fastq.Sequence) error {
	// 使用原子操作更新统计，减少锁竞争
	written := atomic.AddInt64(&bfw.written, 1)

	// 只在需要刷新时获取锁
	if written%100000 == 0 {
		bfw.mu.Lock()
		defer bfw.mu.Unlock()

		// 直接使用fastq.Writer，它已经处理了序列化
		if _, err := bfw.writer.Write(seq); err != nil {
			return err
		}

		return bfw.flush()
	}

	// 大部分情况下不需要锁，直接写入
	if _, err := bfw.writer.Write(seq); err != nil {
		return err
	}

	return nil
}

func (bfw *BufferedFastqWriter) flush() error {
	// 强制同步到磁盘
	return bfw.file.Sync()
}

func (bfw *BufferedFastqWriter) Flush() error {
	bfw.mu.Lock()
	defer bfw.mu.Unlock()
	return bfw.flush()
}

func (bfw *BufferedFastqWriter) Close() error {
	if err := bfw.Flush(); err != nil {
		return err
	}
	return bfw.file.Close()
}

func usage() {
	fmt.Printf("\nProgram: ml2e - FastQ read ID modifier (changes @ML to @E) (use 6-8 threads)\n")
	fmt.Printf("Usage: ml2e [options]\n\n")
	fmt.Printf("Options:\n")
	fmt.Printf("  -fq1        Input fastq file 1 (supports .gz format)\n")
	fmt.Printf("  -fq2        Input fastq file 2 (supports .gz format)\n")
	fmt.Printf("  -outdir     Output directory (required)\n")
	fmt.Printf("  -workers    Number of worker threads (default: 4, recommend 6-8)\n")
	fmt.Printf("  -pigz       Path to pigz executable for compression\n")
	fmt.Printf("  -pigz-decompress  Use pigz for decompression (faster for gzip files)\n")
	fmt.Printf("  -version    Show version information\n\n")
	fmt.Printf("Examples:\n")
	fmt.Printf("  ml2e -fq1 input1.fastq.gz -fq2 input2.fastq.gz -outdir output\n")
	fmt.Printf("  ml2e -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output\n")
	fmt.Printf("  ml2e -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output -pigz /usr/bin/pigz\n")
	fmt.Printf("  ml2e -fq1 input1.fastq.gz -fq2 input2.fastq.gz -outdir output -pigz /usr/bin/pigz -pigz-decompress\n")
	os.Exit(1)
}

// 读取对
type ReadPair struct {
	Seq1 fastq.Sequence
	Seq2 fastq.Sequence
}

// 高性能批量读取流水线模式 - 动态缓存机制
func mainPipelineMode(fq1, fq2, outdir string, numWorkers int, pigzPath string) error {
	// 创建输出目录
	if err := os.MkdirAll(outdir, 0755); err != nil {
		return fmt.Errorf("failed to create output directory: %v", err)
	}

	// 创建输入文件reader
	reader1, file1, err := createReader(fq1)
	if err != nil {
		return err
	}
	defer file1.Close()

	reader2, file2, err := createReader(fq2)
	if err != nil {
		return err
	}
	defer file2.Close()

	// 创建fastq scanner
	scanner1 := fastq.NewScanner(fastq.NewReader(reader1))
	scanner2 := fastq.NewScanner(fastq.NewReader(reader2))

	// 创建输出文件
	basename1 := getBasename(fq1)
	basename2 := getBasename(fq2)
	outFile1Path := filepath.Join(outdir, basename1)
	outFile2Path := filepath.Join(outdir, basename2)

	writer1, outFile1, err := createWriter(outFile1Path)
	if err != nil {
		return err
	}
	defer outFile1.Close()

	writer2, outFile2, err := createWriter(outFile2Path)
	if err != nil {
		return err
	}
	defer outFile2.Close()

	// 创建fastq writer
	w1 := fastq.NewWriter(writer1)
	w2 := fastq.NewWriter(writer2)

	// 创建高性能缓冲写入器，增加缓冲区大小
	bufferedW1 := NewBufferedFastqWriter(w1, outFile1, 64*1024*1024) // 64MB缓冲区（增加100%）
	bufferedW2 := NewBufferedFastqWriter(w2, outFile2, 64*1024*1024) // 64MB缓冲区（增加100%）
	defer bufferedW1.Close()
	defer bufferedW2.Close()

	// 统计变量
	var totalReads int64
	var totalOutputReads int64
	var mu sync.Mutex

	// 动态缓存配置 - 针对gzip文件优化
	const workerBatchSize = 2000 // 每个worker一次处理2000条reads
	const readBatchSize = 50000  // 读取批次大小50K（gzip文件不宜过大）

	// 智能缓存容量计算 - 针对高worker数量优化
	// 写入缓存：worker数量 * 批处理大小 * 缓冲倍数(2-3倍)
	writeCacheMultiplier := 2
	if numWorkers <= 4 {
		writeCacheMultiplier = 3 // 少worker时增加缓冲
	} else if numWorkers >= 8 {
		writeCacheMultiplier = 2 // 多worker时减少缓冲，避免内存过度使用
	}
	writeCacheCapacity := int64(numWorkers) * workerBatchSize * int64(writeCacheMultiplier)

	// 读取队列：worker数量 * 批处理大小 * 缓冲倍数(1-2倍)
	readQueueMultiplier := 1
	if numWorkers <= 4 {
		readQueueMultiplier = 2 // 少worker时增加缓冲
	}
	readQueueCapacity := int64(numWorkers) * workerBatchSize * int64(readQueueMultiplier)

	// 动态写入批次大小 - 基于worker数量优化，减少内存积累
	writeBatchSize := numWorkers * 500 // 每个worker贡献500条reads（减少50%）
	if writeBatchSize < 3000 {
		writeBatchSize = 3000 // 最小保证3000条
	} else if writeBatchSize > 10000 {
		writeBatchSize = 10000 // 最大限制10000条
	}

	// 批量读取配置 - 针对gzip文件优化
	const batchQueueSize = 15 // 适中的队列大小，避免gzip读取阻塞

	// 输出配置信息
	fmt.Printf("Configuration: workers=%d, read_batch=%d, worker_batch=%d, write_batch=%d\n",
		numWorkers, readBatchSize, workerBatchSize, writeBatchSize)

	// 检测gzip文件并输出优化提示
	if strings.HasSuffix(fq1, ".gz") || strings.HasSuffix(fq2, ".gz") {
		fmt.Printf("Gzip files detected - using optimized gzip reading strategy\n")
	}

	// 创建写入缓存池，根据worker批量处理调整容量
	writeCache := make(chan *ProcessResult, writeCacheCapacity)

	// 创建读取队列，根据worker批量处理调整容量
	readQueue := make(chan []ReadPair, readQueueCapacity)

	// 创建控制通道
	stopReading := make(chan bool)
	stopProcessing := make(chan bool)
	stopWriting := make(chan bool)

	// 批量读取配置
	r1BatchQueue := make(chan []fastq.Sequence, batchQueueSize)
	r2BatchQueue := make(chan []fastq.Sequence, batchQueueSize)

	// 启动R1批量读取线程
	r1Done := make(chan bool)
	go func() {
		defer func() {
			r1Done <- true
		}()

		batch := make([]fastq.Sequence, 0, readBatchSize)
		for scanner1.Next() {
			batch = append(batch, scanner1.Seq())
			if len(batch) >= readBatchSize {
				select {
				case r1BatchQueue <- batch:
					batch = make([]fastq.Sequence, 0, readBatchSize)
				case <-stopReading:
					return
				}
			}
		}
		// 发送最后一批
		if len(batch) > 0 {
			select {
			case r1BatchQueue <- batch:
			case <-stopReading:
				return
			}
		}
		close(r1BatchQueue)
	}()

	// 启动R2批量读取线程
	r2Done := make(chan bool)
	go func() {
		defer func() {
			r2Done <- true
		}()

		batch := make([]fastq.Sequence, 0, readBatchSize)
		for scanner2.Next() {
			batch = append(batch, scanner2.Seq())
			if len(batch) >= readBatchSize {
				select {
				case r2BatchQueue <- batch:
					batch = make([]fastq.Sequence, 0, readBatchSize)
				case <-stopReading:
					return
				}
			}
		}
		// 发送最后一批
		if len(batch) > 0 {
			select {
			case r2BatchQueue <- batch:
			case <-stopReading:
				return
			}
		}
		close(r2BatchQueue)
	}()

	// 启动批量配对线程 - 现在发送ReadPair批次而不是单个
	pairDone := make(chan bool)
	go func() {
		defer func() {
			pairDone <- true
		}()

		for {
			select {
			case r1Batch, ok1 := <-r1BatchQueue:
				if !ok1 {
					return
				}
				r2Batch, ok2 := <-r2BatchQueue
				if !ok2 {
					return
				}

				// 批量配对 - 确保R1和R2完全配对
				minLen := len(r1Batch)
				if len(r2Batch) < minLen {
					minLen = len(r2Batch)
				}

				// 将配对的reads分批发送给worker
				readPairBatch := make([]ReadPair, 0, workerBatchSize)
				for i := 0; i < minLen; i++ {
					readPairBatch = append(readPairBatch, ReadPair{Seq1: r1Batch[i], Seq2: r2Batch[i]})

					// 当批次满时发送
					if len(readPairBatch) >= workerBatchSize {
						select {
						case readQueue <- readPairBatch:
							currentReads := atomic.AddInt64(&totalReads, int64(len(readPairBatch)))
							// 每读取2M reads输出进度
							if currentReads%2000000 == 0 {
								fmt.Printf("Read %d reads (batch mode, worker batch size: %d)\n", currentReads, workerBatchSize)
							}
							readPairBatch = make([]ReadPair, 0, workerBatchSize)
						case <-stopReading:
							return
						}
					}
				}

				// 发送最后一批（可能不足workerBatchSize）
				if len(readPairBatch) > 0 {
					select {
					case readQueue <- readPairBatch:
						currentReads := atomic.AddInt64(&totalReads, int64(len(readPairBatch)))
						if currentReads%2000000 == 0 {
							fmt.Printf("Read %d reads (batch mode, worker batch size: %d)\n", currentReads, workerBatchSize)
						}
					case <-stopReading:
						return
					}
				}

				// 检查并报告批次大小差异
				if len(r1Batch) != len(r2Batch) {
					fmt.Printf("Warning: R1 batch size (%d) != R2 batch size (%d), paired %d reads\n",
						len(r1Batch), len(r2Batch), minLen)
				}
			case <-stopReading:
				return
			}
		}
	}()

	// 启动worker线程池 - 现在批量处理
	workerWg := sync.WaitGroup{}
	for i := 0; i < numWorkers; i++ {
		workerWg.Add(1)
		go func(workerID int) {
			defer func() {
				workerWg.Done()
			}()

			for readPairBatch := range readQueue {
				// 批量处理reads
				results := processReadPairBatch(readPairBatch)

				// 批量发送结果 - 所有reads都保留
				for _, result := range results {
					select {
					case writeCache <- result:
					case <-stopProcessing:
						return
					}
				}
			}
		}(i)
	}

	// 启动单函数双协程输出 - 修复通道关闭时机问题
	writeDone := make(chan bool)
	go func() {
		defer func() {
			writeDone <- true
		}()

		// 共享的写入缓存通道
		writeCacheForR1 := make(chan *ProcessResult, writeCacheCapacity)
		writeCacheForR2 := make(chan *ProcessResult, writeCacheCapacity)

		// 使用WaitGroup确保所有协程完成后再关闭通道
		var wg sync.WaitGroup
		wg.Add(2) // R1和R2协程

		// 启动数据分发协程
		go func() {
			defer func() {
				// 等待R1和R2协程完成后再关闭通道
				wg.Wait()
				close(writeCacheForR1)
				close(writeCacheForR2)
			}()

			for {
				select {
				case result := <-writeCache:
					if result == nil {
						// 结束信号
						writeCacheForR1 <- nil
						writeCacheForR2 <- nil
						return
					}

					// 分发R1和R2数据到各自的通道
					select {
					case writeCacheForR1 <- result:
					case <-stopWriting:
						return
					}
					select {
					case writeCacheForR2 <- result:
					case <-stopWriting:
						return
					}
				case <-stopWriting:
					return
				}
			}
		}()

		// 启动R1写入协程
		go func() {
			defer wg.Done()

			r1Batch := make([]fastq.Sequence, 0, writeBatchSize)
			for result := range writeCacheForR1 {
				if result == nil {
					// 结束信号，写入最后一批
					if len(r1Batch) > 0 {
						writeR1Batch(r1Batch, bufferedW1)

						// 更新统计
						mu.Lock()
						totalOutputReads += int64(len(r1Batch))
						mu.Unlock()
					}
					return
				}

				r1Batch = append(r1Batch, result.Seq1)
				if len(r1Batch) >= writeBatchSize {
					writeR1Batch(r1Batch, bufferedW1)

					// 更新统计
					mu.Lock()
					totalOutputReads += int64(len(r1Batch))
					mu.Unlock()

					r1Batch = make([]fastq.Sequence, 0, writeBatchSize)
				}
			}
		}()

		// 启动R2写入协程
		go func() {
			defer wg.Done()

			r2Batch := make([]fastq.Sequence, 0, writeBatchSize)
			for result := range writeCacheForR2 {
				if result == nil {
					// 结束信号，写入最后一批
					if len(r2Batch) > 0 {
						writeR2Batch(r2Batch, bufferedW2)
					}
					return
				}

				r2Batch = append(r2Batch, result.Seq2)
				if len(r2Batch) >= writeBatchSize {
					writeR2Batch(r2Batch, bufferedW2)
					r2Batch = make([]fastq.Sequence, 0, writeBatchSize)
				}
			}
		}()
	}()

	// 等待R1、R2读取线程和配对线程完成
	<-r1Done
	<-r2Done
	<-pairDone

	// 关闭读取相关通道
	close(stopReading)
	close(readQueue)

	// 等待所有worker完成
	workerWg.Wait()
	close(stopProcessing)

	// 发送结束信号给输出线程
	writeCache <- nil

	// 等待输出线程完成
	<-writeDone
	close(stopWriting)

	// 关闭通道
	close(writeCache)

	// 强制刷新缓冲区
	fmt.Println("Flushing all buffers...")
	if err := bufferedW1.Flush(); err != nil {
		fmt.Printf("Error flushing bufferedW1: %v\n", err)
	}
	if err := bufferedW2.Flush(); err != nil {
		fmt.Printf("Error flushing bufferedW2: %v\n", err)
	}

	// 输出统计信息
	statsFile := filepath.Join(outdir, basename2+".stats.txt")
	statsContent := fmt.Sprintf("Total reads processed: %d\nTotal output reads: %d\n",
		totalReads, totalOutputReads)

	if err := os.WriteFile(statsFile, []byte(statsContent), 0644); err != nil {
		return fmt.Errorf("failed to write stats file: %v", err)
	}

	// 输出最终结果
	fmt.Printf("\nRead ID modification pipeline completed. Results saved to: %s\n", statsFile)
	fmt.Printf("Modified files saved to: %s and %s\n", outFile1Path, outFile2Path)
	fmt.Printf("Total reads processed: %d, Total output reads: %d\n",
		totalReads, totalOutputReads)

	// 如果指定了pigz路径，进行压缩
	if pigzPath != "" {
		if err := compressWithPigz(outFile1Path, outFile2Path, pigzPath, numWorkers); err != nil {
			return fmt.Errorf("failed to compress output files: %v", err)
		}
	}

	return nil
}

// 写入R1批次
func writeR1Batch(batch []fastq.Sequence, writer *BufferedFastqWriter) error {
	for _, seq := range batch {
		if err := writer.Write(seq); err != nil {
			log.Printf("Error writing R1: %v", err)
			return err
		}
	}
	return nil
}

// 写入R2批次
func writeR2Batch(batch []fastq.Sequence, writer *BufferedFastqWriter) error {
	for _, seq := range batch {
		if err := writer.Write(seq); err != nil {
			log.Printf("Error writing R2: %v", err)
			return err
		}
	}
	return nil
}

func main() {
	// 定义命令行参数
	var fq1, fq2, outdir, pigzPath string
	var numWorkers int
	var showVersion bool
	var usePigzDecompress bool

	flag.StringVar(&fq1, "fq1", "", "Input fastq file 1 (supports .gz format)")
	flag.StringVar(&fq2, "fq2", "", "Input fastq file 2 (supports .gz format)")
	flag.StringVar(&outdir, "outdir", "", "Output directory")
	flag.IntVar(&numWorkers, "workers", 4, "Number of worker threads")
	flag.StringVar(&pigzPath, "pigz", "", "Path to pigz executable for compression")
	flag.BoolVar(&showVersion, "version", false, "Show version information")
	flag.BoolVar(&usePigzDecompress, "pigz-decompress", false, "Use pigz for decompression (faster for gzip files)")

	// 解析命令行参数
	flag.Parse()

	// 显示版本信息
	if showVersion {
		fmt.Printf("ml2e v%s (Build: %s)\n", Version, BuildTime)
		os.Exit(0)
	}

	// 验证必需参数
	if fq1 == "" || fq2 == "" || outdir == "" {
		fmt.Println("Error: -fq1, -fq2 and -outdir are required")
		usage()
	}

	if numWorkers < 1 {
		fmt.Println("Error: workers must be at least 1")
		usage()
	}

	// 执行高性能批量读取流水线模式
	fmt.Println("Starting read ID modification pipeline...")

	// 如果启用pigz解压缩且指定了pigz路径
	if usePigzDecompress && pigzPath != "" {
		fmt.Println("Using pigz for decompression...")

		// 解压缩gzip文件
		if strings.HasSuffix(fq1, ".gz") {
			decompressedFq1, err := decompressWithPigz(fq1, pigzPath, numWorkers)
			if err != nil {
				log.Fatalf("Failed to decompress %s: %v", fq1, err)
			}
			fq1 = decompressedFq1
			defer os.Remove(decompressedFq1) // 处理完成后删除临时文件
		}

		if strings.HasSuffix(fq2, ".gz") {
			decompressedFq2, err := decompressWithPigz(fq2, pigzPath, numWorkers)
			if err != nil {
				log.Fatalf("Failed to decompress %s: %v", fq2, err)
			}
			fq2 = decompressedFq2
			defer os.Remove(decompressedFq2) // 处理完成后删除临时文件
		}
	}

	err := mainPipelineMode(fq1, fq2, outdir, numWorkers, pigzPath)
	if err != nil {
		log.Fatalf("Error: %v", err)
	}
}
