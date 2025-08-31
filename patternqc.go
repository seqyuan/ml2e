package main

import (
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/seqyuan/annogene/io/fastq"
)

const (
	Version   = "1.0.0"
	BuildTime = "2024-08-30"
)

// 处理结果结构
type ProcessResult struct {
	Index      int // 读取顺序索引
	Seq1       fastq.Sequence
	Seq2       fastq.Sequence
	HasPattern bool
	Keep       bool
}

// 线程安全的顺序队列 - 高性能版本
type OrderedQueue struct {
	mu           sync.Mutex
	results      map[int]*ProcessResult
	nextIndex    int
	cond         *sync.Cond
	closed       bool
	maxQueueSize int // 新增：限制队列最大大小
}

func NewOrderedQueue() *OrderedQueue {
	q := &OrderedQueue{
		results:      make(map[int]*ProcessResult),
		nextIndex:    0,
		closed:       false,
		maxQueueSize: 5000000, // 增加到500万个结果，大幅提高并发性能
	}
	q.cond = sync.NewCond(&q.mu)
	return q
}

// 添加结果到队列 - 无阻塞版本
func (q *OrderedQueue) Add(result *ProcessResult) {
	q.mu.Lock()
	defer q.mu.Unlock()

	if q.closed {
		return
	}

	// 移除队列大小限制，让worker完全不被阻塞
	q.results[result.Index] = result
	q.cond.Signal()
}

// 获取下一个可输出的结果 - 优化版本
func (q *OrderedQueue) GetNext() *ProcessResult {
	q.mu.Lock()
	defer q.mu.Unlock()

	for {
		if result, exists := q.results[q.nextIndex]; exists {
			delete(q.results, q.nextIndex)
			q.nextIndex++
			return result
		}

		// 检查是否所有工作都完成了
		if q.closed && len(q.results) == 0 {
			return nil // 表示结束
		}

		// 使用更短的超时时间，减少阻塞
		done := make(chan struct{})
		go func() {
			q.cond.Wait()
			close(done)
		}()

		select {
		case <-done:
			// 继续循环
		case <-time.After(100 * time.Microsecond): // 减少到100微秒
			// 超时，检查是否有数据丢失
			if len(q.results) > 0 {
				// 找到最小的index
				minIndex := q.nextIndex
				for idx := range q.results {
					if idx < minIndex {
						minIndex = idx
					}
				}
				if minIndex < q.nextIndex {
					// 有数据丢失，调整nextIndex
					fmt.Printf("Warning: Data gap detected, adjusting nextIndex from %d to %d\n", q.nextIndex, minIndex)
					q.nextIndex = minIndex
				}
			}
		}
	}
}

// 关闭队列
func (q *OrderedQueue) Close() {
	q.mu.Lock()
	defer q.mu.Unlock()

	q.closed = true
	q.cond.Broadcast()
}

// 检查队列状态
func (q *OrderedQueue) IsEmpty() bool {
	q.mu.Lock()
	defer q.mu.Unlock()
	return len(q.results) == 0
}

// 获取队列长度
func (q *OrderedQueue) Len() int {
	q.mu.Lock()
	defer q.mu.Unlock()
	return len(q.results)
}

// 获取队列详细信息
func (q *OrderedQueue) GetStatus() (int, int, int) {
	q.mu.Lock()
	defer q.mu.Unlock()

	minIndex := q.nextIndex
	maxIndex := q.nextIndex
	for idx := range q.results {
		if idx < minIndex {
			minIndex = idx
		}
		if idx > maxIndex {
			maxIndex = idx
		}
	}

	return len(q.results), q.nextIndex, maxIndex
}

// 强制清理队列中的旧数据 - 优化版本
func (q *OrderedQueue) Cleanup() {
	q.mu.Lock()
	defer q.mu.Unlock()

	// 只在队列过大时才清理
	if len(q.results) > q.maxQueueSize*95/100 { // 超过95%时才清理
		// 找到需要清理的索引范围
		toDelete := make([]int, 0)
		for idx := range q.results {
			if idx < q.nextIndex-q.maxQueueSize/20 { // 只清理最旧的5%
				toDelete = append(toDelete, idx)
			}
		}

		// 删除旧数据
		for _, idx := range toDelete {
			delete(q.results, idx)
		}

		if len(toDelete) > 0 {
			fmt.Printf("Cleaned up %d old results from queue (queue size: %d)\n", len(toDelete), len(q.results))
		}
	}
}

// 反向互补序列
func reverseComplement(pattern string) string {
	complement := map[byte]byte{
		'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
		'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
		'N': 'N', 'n': 'n',
	}

	result := make([]byte, len(pattern))
	for i := len(pattern) - 1; i >= 0; i-- {
		if comp, exists := complement[pattern[i]]; exists {
			result[len(pattern)-1-i] = comp
		} else {
			result[len(pattern)-1-i] = pattern[i]
		}
	}
	return string(result)
}

// 检查序列是否包含pattern或其反向互补序列
func containsPattern(sequence, pattern string) bool {
	revComp := reverseComplement(pattern)
	return strings.Contains(sequence, pattern) || strings.Contains(sequence, revComp)
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
		// 使用gzip reader
		gz, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, nil, fmt.Errorf("failed to create gzip reader for %s: %v", filePath, err)
		}
		reader = gz
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

// 处理单个read pair的worker函数
func processReadPair(seq1, seq2 fastq.Sequence, index int, pattern string, keepEveryN int, patternCount *int, mu *sync.Mutex) *ProcessResult {
	hasPattern := containsPattern(string(seq2.Letters), pattern)

	var keep bool
	if hasPattern {
		mu.Lock()
		*patternCount++
		currentCount := *patternCount
		mu.Unlock()

		keep = currentCount%keepEveryN == 0
	} else {
		keep = true // 不包含pattern的read直接保留
	}

	return &ProcessResult{
		Index:      index,
		Seq1:       seq1,
		Seq2:       seq2,
		HasPattern: hasPattern,
		Keep:       keep,
	}
}

// Filter模式：过滤fq1和fq2，保留指定比例的pattern reads
func filterMode(fq1, fq2, pattern, outdir string, percent int, numWorkers int, pigzPath string) error {
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

	// 准备输出文件路径 - 保持原始扩展名，但移除.gz
	basename1 := getBasename(fq1)
	basename2 := getBasename(fq2)

	// 输出文件使用.fq扩展名，不压缩
	output1 := filepath.Join(outdir, basename1)
	output2 := filepath.Join(outdir, basename2)

	// 创建输出文件
	writer1, outFile1, err := createWriter(output1)
	if err != nil {
		return err
	}
	defer outFile1.Close()

	writer2, outFile2, err := createWriter(output2)
	if err != nil {
		return err
	}
	defer outFile2.Close()

	// 创建fastq writer
	w1 := fastq.NewWriter(writer1)
	w2 := fastq.NewWriter(writer2)

	// 创建顺序队列和同步原语
	queue := NewOrderedQueue()
	var wg sync.WaitGroup
	var mu sync.Mutex

	totalReads := 0
	patternReads := 0
	keptPatternReads := 0
	totalOutputReads := 0 // 新增：统计所有输出的reads

	// 计算保留的pattern reads数量
	keepEveryN := 100 / percent

	// 创建worker池 - 无限制缓冲区
	workChan := make(chan struct {
		seq1  fastq.Sequence
		seq2  fastq.Sequence
		index int
	}, numWorkers*1000) // 极大增加缓冲区，确保worker完全不被阻塞

	// 启动worker goroutines - 高性能版本
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func(workerID int) {
			defer wg.Done()
			processedCount := 0
			for work := range workChan {
				result := processReadPair(work.seq1, work.seq2, work.index, pattern, keepEveryN, &patternReads, &mu)
				queue.Add(result)
				processedCount++
				
				// 每处理10000个reads输出一次worker状态
				if processedCount%10000 == 0 {
					fmt.Printf("Worker %d processed %d reads\n", workerID, processedCount)
				}
			}
			fmt.Printf("Worker %d finished, total processed: %d\n", workerID, processedCount)
		}(i)
	}

	// 启动输出goroutine - 高性能版本
	outputDone := make(chan bool)
	go func() {
		defer func() { outputDone <- true }()

		// 增加批量写入大小，提高I/O效率
		batchSize := 50000 // 从5000大幅增加到50000
		batch := make([]*ProcessResult, 0, batchSize)
		cleanupCounter := 0

		for {
			result := queue.GetNext()
			if result == nil {
				// 队列已关闭，处理剩余批次
				if len(batch) > 0 {
					writeBatch(batch, w1, w2, &keptPatternReads, &totalOutputReads, &mu)
				}
				break
			}

			batch = append(batch, result)

			// 当批次满了时写入
			if len(batch) >= batchSize {
				writeBatch(batch, w1, w2, &keptPatternReads, &totalOutputReads, &mu)
				batch = batch[:0] // 清空批次

				// 大幅减少清理频率
				cleanupCounter++
				if cleanupCounter >= 1000 { // 每1000个批次清理一次
					queue.Cleanup()
					cleanupCounter = 0
				}
			}
		}
	}()

	// 读取并分发工作 - 高性能版本
	index := 0
	
	for scanner1.Next() && scanner2.Next() {
		seq1 := scanner1.Seq()
		seq2 := scanner2.Seq()
		totalReads++

		// 发送工作到worker池 - 无阻塞版本
		workChan <- struct {
			seq1  fastq.Sequence
			seq2  fastq.Sequence
			index int
		}{seq1, seq2, index}

		index++

		// 每处理100000个reads输出一次进度，包括队列状态和内存使用
		if totalReads%100000 == 0 {
			queueLen, nextIdx, maxIdx := queue.GetStatus()

			// 获取内存使用情况
			var m runtime.MemStats
			runtime.ReadMemStats(&m)

			fmt.Printf("Processed %d reads, Queue: %d, Next: %d, Max: %d, Memory: %.2f GB\n",
				totalReads, queueLen, nextIdx, maxIdx, float64(m.Alloc)/1024/1024/1024)

			// 如果内存使用过高，强制垃圾回收（调整到48GB）
			if m.Alloc > 48*1024*1024*1024 { // 超过48GB
				fmt.Printf("High memory usage detected (%.2f GB), forcing garbage collection...\n",
					float64(m.Alloc)/1024/1024/1024)
				runtime.GC()
				queue.Cleanup()
			}
		}
	}

	// 关闭工作通道
	close(workChan)

	// 等待所有worker完成
	wg.Wait()

	// 关闭队列
	queue.Close()

	// 等待输出完成，根据数据量动态调整超时时间
	timeout := time.Duration(totalReads/1000000+1) * 60 * time.Second // 每百万reads增加1分钟
	if timeout < 30*time.Second {
		timeout = 30 * time.Second
	}
	if timeout > 600*time.Second {
		timeout = 600 * time.Second // 最大10分钟
	}

	fmt.Printf("Waiting for output completion (timeout: %v)...\n", timeout)
	select {
	case <-outputDone:
		fmt.Println("Output completed successfully")
	case <-time.After(timeout):
		fmt.Printf("Warning: Output timeout after %v, but processing completed\n", timeout)
		queueLen, nextIdx, maxIdx := queue.GetStatus()
		fmt.Printf("Final queue status - Length: %d, Next: %d, Max: %d\n", queueLen, nextIdx, maxIdx)
	}

	if err := scanner1.Error(); err != nil {
		return fmt.Errorf("error reading fq1: %v", err)
	}
	if err := scanner2.Error(); err != nil {
		return fmt.Errorf("error reading fq2: %v", err)
	}

	// 计算最终比例
	finalRatio := 0.0
	if totalReads > 0 {
		finalRatio = float64(keptPatternReads) / float64(totalReads) * 100.0
	}

	// 输出统计结果
	outputFile := filepath.Join(outdir, fmt.Sprintf("%s.stats.txt", basename2))
	output, err := os.Create(outputFile)
	if err != nil {
		return fmt.Errorf("failed to create ratio file: %v", err)
	}
	defer output.Close()

	fmt.Fprintf(output, "Total reads: %d\n", totalReads)
	fmt.Fprintf(output, "Original pattern reads: %d\n", patternReads)
	fmt.Fprintf(output, "Kept pattern reads: %d\n", keptPatternReads)
	fmt.Fprintf(output, "Total output reads: %d\n", totalOutputReads) // 新增：总输出reads数
	fmt.Fprintf(output, "Final ratio: %.2f%%\n", finalRatio)

	fmt.Printf("Filter mode completed. Results saved to: %s\n", outputFile)
	fmt.Printf("Filtered files saved to: %s and %s\n", output1, output2)
	fmt.Printf("Total reads: %d, Kept pattern reads: %d, Total output reads: %d, Final ratio: %.2f%%\n",
		totalReads, keptPatternReads, totalOutputReads, finalRatio)

	// 如果指定了pigz路径，则压缩输出文件
	if pigzPath != "" {
		if err := compressWithPigz(output1, output2, pigzPath, numWorkers); err != nil {
			log.Printf("Warning: Failed to compress output files with pigz: %v", err)
		} else {
			fmt.Printf("Output files compressed with pigz using %d threads\n", numWorkers)
		}
	}

	return nil
}

// 使用pigz压缩文件
func compressWithPigz(file1, file2, pigzPath string, numWorkers int) error {
	// 检查pigz路径是否存在
	if _, err := os.Stat(pigzPath); os.IsNotExist(err) {
		return fmt.Errorf("pigz not found at path: %s", pigzPath)
	}

	// 压缩第一个文件
	if err := compressFileWithPigz(file1, pigzPath, numWorkers); err != nil {
		return fmt.Errorf("failed to compress %s: %v", file1, err)
	}

	// 压缩第二个文件
	if err := compressFileWithPigz(file2, pigzPath, numWorkers); err != nil {
		return fmt.Errorf("failed to compress %s: %v", file2, err)
	}

	return nil
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

// 批量写入函数
func writeBatch(batch []*ProcessResult, w1, w2 fastq.Writer, keptPatternReads *int, totalOutputReads *int, mu *sync.Mutex) {
	for _, result := range batch {
		if result.Keep {
			// 写入fq1 - 确保写入完整的FASTQ记录（包括质量值）
			if _, err := w1.Write(result.Seq1); err != nil {
				log.Printf("error writing to fq1: %v", err)
				return
			}
			// 写入fq2 - 确保写入完整的FASTQ记录（包括质量值）
			if _, err := w2.Write(result.Seq2); err != nil {
				log.Printf("error writing to fq2: %v", err)
				return
			}

			// 统计所有输出的reads
			mu.Lock()
			*totalOutputReads++
			if result.HasPattern {
				*keptPatternReads++
			}
			mu.Unlock()
		}
	}
}

func usage() {
	fmt.Printf("\nProgram: patternqc - FastQ filtering based on TSO pattern (Multi-threaded, Optimized)\n")
	fmt.Printf("Usage: patternqc [options]\n\n")
	fmt.Printf("Options:\n")
	fmt.Printf("  -fq1        Input fastq file 1 (supports .gz format)\n")
	fmt.Printf("  -fq2        Input fastq file 2 (supports .gz format)\n")
	fmt.Printf("  -pattern    Pattern to search for (default: AGCAGTGGTATCAACGCAGAGTACA)\n")
	fmt.Printf("  -outdir     Output directory (required)\n")
	fmt.Printf("  -percent    Percentage of pattern reads to keep (0-100, default: 5)\n")
	fmt.Printf("  -workers    Number of worker threads (default: 4)\n")
	fmt.Printf("  -pigz       Path to pigz executable for compression\n")
	fmt.Printf("  -version    Show version information\n\n")
	fmt.Printf("Examples:\n")
	fmt.Printf("  patternqc -fq1 input1.fastq.gz -fq2 input2.fastq.gz -outdir output\n")
	fmt.Printf("  patternqc -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output -percent 10 -workers 8\n")
	fmt.Printf("  patternqc -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output -pigz /usr/bin/pigz\n")
	os.Exit(1)
}

func main() {
	// 定义命令行参数
	var fq1, fq2, pattern, outdir, pigzPath string
	var percent, numWorkers int
	var showVersion bool

	flag.StringVar(&fq1, "fq1", "", "Input fastq file 1 (supports .gz format)")
	flag.StringVar(&fq2, "fq2", "", "Input fastq file 2 (supports .gz format)")
	flag.StringVar(&pattern, "pattern", "AGCAGTGGTATCAACGCAGAGTACA", "Pattern to search for")
	flag.StringVar(&outdir, "outdir", "", "Output directory")
	flag.IntVar(&percent, "percent", 5, "Percentage of pattern reads to keep (0-100)")
	flag.IntVar(&numWorkers, "workers", 4, "Number of worker threads")
	flag.StringVar(&pigzPath, "pigz", "", "Path to pigz executable for compression")
	flag.BoolVar(&showVersion, "version", false, "Show version information")

	// 解析命令行参数
	flag.Parse()

	// 显示版本信息
	if showVersion {
		fmt.Printf("PatternQC v%s (Build: %s)\n", Version, BuildTime)
		fmt.Printf("GitHub: https://github.com/seqyuan/patternqc\n")
		os.Exit(0)
	}

	// 验证必需参数
	if fq1 == "" || fq2 == "" || outdir == "" {
		fmt.Println("Error: -fq1, -fq2 and -outdir are required")
		usage()
	}

	if percent < 0 || percent > 100 {
		fmt.Println("Error: percent must be between 0 and 100")
		usage()
	}

	if numWorkers < 1 {
		fmt.Println("Error: workers must be at least 1")
		usage()
	}

	// 执行过滤模式
	err := filterMode(fq1, fq2, pattern, outdir, percent, numWorkers, pigzPath)
	if err != nil {
		log.Fatalf("Error: %v", err)
	}
}
