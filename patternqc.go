package main

import (
	"bytes"
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
	"sync/atomic"
	"time"

	"github.com/edsrzf/mmap-go"
	"github.com/seqyuan/annogene/io/fastq"
)

const (
	Version   = "5.1.0"
	BuildTime = "2024-12-19"
)

// 简化的结果结构 - 只需要配对信息
type ProcessResult struct {
	Seq1       fastq.Sequence
	Seq2       fastq.Sequence
	HasPattern bool
	Keep       bool
}

// 简化的队列 - 不需要顺序控制
type SimpleQueue struct {
	mu      sync.Mutex
	results []*ProcessResult
	closed  bool
}

func NewSimpleQueue() *SimpleQueue {
	return &SimpleQueue{
		results: make([]*ProcessResult, 0, 100000), // 预分配10万容量
		closed:  false,
	}
}

// 添加结果到队列 - 无顺序要求
func (q *SimpleQueue) Add(result *ProcessResult) {
	q.mu.Lock()
	defer q.mu.Unlock()

	if q.closed {
		return
	}

	q.results = append(q.results, result)
}

// 获取所有结果并清空队列
func (q *SimpleQueue) GetAll() []*ProcessResult {
	q.mu.Lock()
	defer q.mu.Unlock()

	if len(q.results) == 0 {
		if q.closed {
			return nil // 队列已关闭且为空
		}
		return []*ProcessResult{} // 返回空切片，表示暂时没有数据
	}

	// 返回所有结果并清空队列
	results := q.results
	q.results = make([]*ProcessResult, 0, 100000)
	return results
}

// 关闭队列
func (q *SimpleQueue) Close() {
	q.mu.Lock()
	defer q.mu.Unlock()
	q.closed = true
}

// 获取队列长度
func (q *SimpleQueue) Len() int {
	q.mu.Lock()
	defer q.mu.Unlock()
	return len(q.results)
}

// 检查队列是否已关闭
func (q *SimpleQueue) IsClosed() bool {
	q.mu.Lock()
	defer q.mu.Unlock()
	return q.closed
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
func processReadPair(seq1, seq2 fastq.Sequence, pattern string, keepEveryN int, patternCount *int, mu *sync.Mutex) *ProcessResult {
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

	// 创建高性能缓冲写入器
	bufferedW1 := NewBufferedFastqWriter(w1, outFile1, 8*1024*1024) // 增加到8MB缓冲区
	bufferedW2 := NewBufferedFastqWriter(w2, outFile2, 8*1024*1024) // 增加到8MB缓冲区
	defer bufferedW1.Close()
	defer bufferedW2.Close()

	// 创建顺序队列和同步原语
	queue := NewSimpleQueue()
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
		seq1 fastq.Sequence
		seq2 fastq.Sequence
	}, numWorkers*2000) // 进一步增加缓冲区，确保worker完全不被阻塞

	// 启动worker goroutines - 高性能版本
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func(workerID int) {
			defer wg.Done()
			processedCount := 0
			for work := range workChan {
				result := processReadPair(work.seq1, work.seq2, pattern, keepEveryN, &patternReads, &mu)
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

	// 启动多个输出goroutine - 真正的多线程输出
	outputDone := make(chan bool)
	numOutputWorkers := 4 // 使用4个输出线程

	for i := 0; i < numOutputWorkers; i++ {
		go func(outputWorkerID int) {
			defer func() {
				if outputWorkerID == 0 {
					outputDone <- true
				}
			}()

			// 每个输出线程处理自己的批次
			batchSize := 1000000 // 恢复到100万个结果
			batch := make([]*ProcessResult, 0, batchSize)

			for {
				results := queue.GetAll()
				if results == nil {
					// 队列已关闭且为空，处理剩余批次
					if len(batch) > 0 {
						writeBatchOptimized(batch, bufferedW1, bufferedW2, &keptPatternReads, &totalOutputReads, &mu)
						fmt.Printf("Output worker %d wrote final batch of %d results\n", outputWorkerID, len(batch))
					}
					break
				}

				// 如果有数据，添加到批次
				if len(results) > 0 {
					batch = append(batch, results...)
				}

				// 当批次满了时写入
				if len(batch) >= batchSize {
					writeBatchOptimized(batch, bufferedW1, bufferedW2, &keptPatternReads, &totalOutputReads, &mu)
					fmt.Printf("Output worker %d wrote batch of %d results\n", outputWorkerID, len(batch))
					batch = batch[:0] // 清空批次
				}

				// 如果队列已关闭且批次不为空，立即写入
				if queue.IsClosed() && len(batch) > 0 {
					writeBatchOptimized(batch, bufferedW1, bufferedW2, &keptPatternReads, &totalOutputReads, &mu)
					fmt.Printf("Output worker %d wrote remaining batch of %d results\n", outputWorkerID, len(batch))
					batch = batch[:0]
				}
			}
		}(i)
	}

	// 读取并分发工作 - 高性能版本
	index := 0

	for scanner1.Next() && scanner2.Next() {
		seq1 := scanner1.Seq()
		seq2 := scanner2.Seq()
		totalReads++

		// 发送工作到worker池 - 无阻塞版本
		workChan <- struct {
			seq1 fastq.Sequence
			seq2 fastq.Sequence
		}{seq1, seq2}

		index++

		// 每处理100000个reads输出一次进度，包括队列状态和内存使用
		if totalReads%100000 == 0 {
			queueLen := queue.Len()

			// 获取内存使用情况
			var m runtime.MemStats
			runtime.ReadMemStats(&m)

			fmt.Printf("Processed %d reads, Queue: %d, Memory: %.2f GB\n",
				totalReads, queueLen, float64(m.Alloc)/1024/1024/1024)

			// 如果内存使用过高，强制垃圾回收（调整到50GB）
			if m.Alloc > 50*1024*1024*1024 { // 超过50GB
				fmt.Printf("High memory usage detected (%.2f GB), forcing garbage collection...\n",
					float64(m.Alloc)/1024/1024/1024)
				runtime.GC()
			}
		}
	}

	// 关闭工作通道
	close(workChan)

	// 等待所有worker完成
	wg.Wait()

	// 关闭队列
	queue.Close()

	// 等待输出完成，确保所有数据都写入
	fmt.Printf("Waiting for output completion...\n")

	// 等待所有输出goroutine完成
	select {
	case <-outputDone:
		fmt.Println("Output completed successfully")
	case <-time.After(30 * time.Second):
		fmt.Printf("Warning: Output timeout after 30s, checking final status...\n")
		queueLen := queue.Len()
		fmt.Printf("Final queue status - Length: %d, Closed: %v\n", queueLen, queue.IsClosed())

		// 强制等待更长时间，确保数据写入
		select {
		case <-outputDone:
			fmt.Println("Output completed after extended wait")
		case <-time.After(60 * time.Second):
			fmt.Printf("Error: Output timeout after 90s total, some data may be lost\n")
		}
	}

	// 强制刷新所有缓冲区
	fmt.Println("Flushing all buffers...")
	if err := bufferedW1.Flush(); err != nil {
		fmt.Printf("Error flushing bufferedW1: %v\n", err)
	}
	if err := bufferedW2.Flush(); err != nil {
		fmt.Printf("Error flushing bufferedW2: %v\n", err)
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

// 分离的写入函数 - 每个文件使用独立线程
func writeBatchSeparate(batch []*ProcessResult, w1, w2 fastq.Writer, keptPatternReads *int, totalOutputReads *int, mu *sync.Mutex) {
	// 分离R1和R2的数据
	var r1Data, r2Data []fastq.Sequence

	for _, result := range batch {
		if result.Keep {
			r1Data = append(r1Data, result.Seq1)
			r2Data = append(r2Data, result.Seq2)
		}
	}

	// 使用goroutine并行写入两个文件
	var wg sync.WaitGroup
	wg.Add(2)

	// 写入R1文件
	go func() {
		defer wg.Done()
		for _, seq := range r1Data {
			if _, err := w1.Write(seq); err != nil {
				log.Printf("error writing to fq1: %v", err)
				return
			}
		}
	}()

	// 写入R2文件
	go func() {
		defer wg.Done()
		for _, seq := range r2Data {
			if _, err := w2.Write(seq); err != nil {
				log.Printf("error writing to fq2: %v", err)
				return
			}
		}
	}()

	// 等待两个写入完成
	wg.Wait()

	// 统计输出
	mu.Lock()
	*totalOutputReads += len(r1Data)
	for _, result := range batch {
		if result.Keep && result.HasPattern {
			*keptPatternReads++
		}
	}
	mu.Unlock()
}

// 高性能缓冲写入器
type BufferedFastqWriter struct {
	mu      sync.Mutex
	buffer  []byte
	writer  fastq.Writer
	file    *os.File
	maxSize int
	written int64
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
	bfw.mu.Lock()
	defer bfw.mu.Unlock()

	// 直接使用fastq.Writer，它已经处理了序列化
	if _, err := bfw.writer.Write(seq); err != nil {
		return err
	}

	bfw.written++

	// 定期刷新缓冲区
	if bfw.written%10000 == 0 {
		return bfw.flush()
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

// 高性能分离写入函数
func writeBatchOptimized(batch []*ProcessResult, w1, w2 *BufferedFastqWriter, keptPatternReads *int, totalOutputReads *int, mu *sync.Mutex) {
	// 分离R1和R2的数据
	var r1Data, r2Data []fastq.Sequence

	for _, result := range batch {
		if result.Keep {
			r1Data = append(r1Data, result.Seq1)
			r2Data = append(r2Data, result.Seq2)
		}
	}

	// 使用goroutine并行写入两个文件
	var wg sync.WaitGroup
	var r1Err, r2Err error
	wg.Add(2)

	// 写入R1文件
	go func() {
		defer wg.Done()
		for _, seq := range r1Data {
			if err := w1.Write(seq); err != nil {
				r1Err = err
				log.Printf("error writing to fq1: %v", err)
				return
			}
		}
	}()

	// 写入R2文件
	go func() {
		defer wg.Done()
		for _, seq := range r2Data {
			if err := w2.Write(seq); err != nil {
				r2Err = err
				log.Printf("error writing to fq2: %v", err)
				return
			}
		}
	}()

	// 等待两个写入完成
	wg.Wait()

	// 检查错误
	if r1Err != nil || r2Err != nil {
		log.Printf("Write errors - R1: %v, R2: %v", r1Err, r2Err)
		return
	}

	// 统计输出
	mu.Lock()
	*totalOutputReads += len(r1Data)
	for _, result := range batch {
		if result.Keep && result.HasPattern {
			*keptPatternReads++
		}
	}
	mu.Unlock()
}

// 借鉴fastp的内存池
type MemoryPool struct {
	pool sync.Pool
	size int
}

func NewMemoryPool(size int) *MemoryPool {
	return &MemoryPool{
		pool: sync.Pool{
			New: func() interface{} {
				return make([]byte, 0, size)
			},
		},
		size: size,
	}
}

func (mp *MemoryPool) Get() []byte {
	return mp.pool.Get().([]byte)
}

func (mp *MemoryPool) Put(buf []byte) {
	// 重置缓冲区
	buf = buf[:0]
	mp.pool.Put(buf)
}

// 借鉴fastp的线程池设计
type ThreadPool struct {
	workers    int
	taskQueue  chan interface{}
	workerFunc func(interface{})
	wg         sync.WaitGroup
}

func NewThreadPool(workers int, workerFunc func(interface{})) *ThreadPool {
	tp := &ThreadPool{
		workers:    workers,
		taskQueue:  make(chan interface{}, workers*100),
		workerFunc: workerFunc,
	}

	for i := 0; i < workers; i++ {
		tp.wg.Add(1)
		go tp.worker()
	}

	return tp
}

func (tp *ThreadPool) worker() {
	defer tp.wg.Done()
	for task := range tp.taskQueue {
		tp.workerFunc(task)
	}
}

func (tp *ThreadPool) Submit(task interface{}) {
	tp.taskQueue <- task
}

func (tp *ThreadPool) Close() {
	close(tp.taskQueue)
	tp.wg.Wait()
}

// 借鉴fastp的智能写入器
type SmartFastqWriter struct {
	mu         sync.Mutex
	buffers    map[string]*bytes.Buffer
	writers    map[string]*BufferedFastqWriter
	batchSize  int
	flushTimer *time.Timer
	memoryPool *MemoryPool
}

func NewSmartFastqWriter(w1, w2 *BufferedFastqWriter, batchSize int) *SmartFastqWriter {
	sfw := &SmartFastqWriter{
		buffers: map[string]*bytes.Buffer{
			"R1": bytes.NewBuffer(make([]byte, 0, batchSize)),
			"R2": bytes.NewBuffer(make([]byte, 0, batchSize)),
		},
		writers: map[string]*BufferedFastqWriter{
			"R1": w1,
			"R2": w2,
		},
		batchSize:  batchSize,
		memoryPool: NewMemoryPool(batchSize),
	}

	// 设置定时刷新
	sfw.flushTimer = time.AfterFunc(100*time.Millisecond, func() {
		sfw.FlushAll()
	})

	return sfw
}

func (sfw *SmartFastqWriter) Write(fileType string, seq fastq.Sequence) error {
	sfw.mu.Lock()
	defer sfw.mu.Unlock()

	buffer := sfw.buffers[fileType]

	// 序列化序列
	serialized := sfw.serializeSequence(seq)
	buffer.Write(serialized)

	// 检查是否需要刷新
	if buffer.Len() >= sfw.batchSize {
		return sfw.flushFile(fileType)
	}

	return nil
}

func (sfw *SmartFastqWriter) flushFile(fileType string) error {
	buffer := sfw.buffers[fileType]
	writer := sfw.writers[fileType]

	if buffer.Len() == 0 {
		return nil
	}

	// 直接写入文件，绕过fastq.Writer
	data := buffer.Bytes()
	_, err := writer.file.Write(data)
	if err != nil {
		return err
	}

	// 重置缓冲区
	buffer.Reset()

	return nil
}

func (sfw *SmartFastqWriter) FlushAll() error {
	sfw.mu.Lock()
	defer sfw.mu.Unlock()

	for fileType := range sfw.buffers {
		if err := sfw.flushFile(fileType); err != nil {
			return err
		}
	}

	// 重置定时器
	sfw.flushTimer.Reset(100 * time.Millisecond)

	return nil
}

func (sfw *SmartFastqWriter) serializeSequence(seq fastq.Sequence) []byte {
	// 使用内存池获取缓冲区
	buf := sfw.memoryPool.Get()
	defer sfw.memoryPool.Put(buf)

	// 序列化逻辑
	// ... 实现序列化

	return buf
}

func usage() {
	fmt.Printf("\nProgram: patternqc - FastQ filtering based on sequence pattern (Multi-threaded, Optimized)\n")
	fmt.Printf("Usage: patternqc [options]\n\n")
	fmt.Printf("Options:\n")
	fmt.Printf("  -fq1        Input fastq file 1 (supports .gz format)\n")
	fmt.Printf("  -fq2        Input fastq file 2 (supports .gz format)\n")
	fmt.Printf("  -pattern    Pattern to search for (default: AGCAGTGGTATCAACGCAGAGTACA)\n")
	fmt.Printf("  -outdir     Output directory (required)\n")
	fmt.Printf("  -percent    Percentage of pattern reads to keep (0-100, default: 5)\n")
	fmt.Printf("  -workers    Number of worker threads (default: 4)\n")
	fmt.Printf("  -mmap       Use memory-mapped file reading for better performance\n")
	fmt.Printf("  -pigz       Path to pigz executable for compression\n")
	fmt.Printf("  -version    Show version information\n\n")
	fmt.Printf("Examples:\n")
	fmt.Printf("  patternqc -fq1 input1.fastq.gz -fq2 input2.fastq.gz -outdir output\n")
	fmt.Printf("  patternqc -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output -percent 10 -workers 8\n")
	fmt.Printf("  patternqc -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output -pigz /usr/bin/pigz\n")
	fmt.Printf("  patternqc -fq1 input1.fastq -fq2 input2.fastq.gz -outdir output -mmap\n")
	os.Exit(1)
}

// 新的流水线架构 - 分离读取、处理、写入
type PipelineArchitecture struct {
	// 配置参数
	batchSize       int // 每批次reads数量 (5M)
	maxCacheBatches int // 最大缓存批次数量 (300)

	// 读取相关
	scanner1, scanner2 *fastq.Scanner

	// 处理相关
	workers     int
	workChan    chan *ReadBatch
	processChan chan *ProcessedBatch

	// 写入相关
	writeChan chan *WriteBatch
	writeDone chan bool

	// 统计
	totalReads       int64
	patternReads     int64
	keptPatternReads int64
	totalOutputReads int64

	// 同步
	wg sync.WaitGroup
	mu sync.Mutex
}

// 读取批次
type ReadBatch struct {
	ID       int
	Reads    []ReadPair
	Complete bool
}

// 处理后的批次
type ProcessedBatch struct {
	ID       int
	Results  []*ProcessResult
	Complete bool
}

// 写入批次
type WriteBatch struct {
	ID       int
	R1Data   []fastq.Sequence
	R2Data   []fastq.Sequence
	Complete bool
}

// 读取对
type ReadPair struct {
	Seq1 fastq.Sequence
	Seq2 fastq.Sequence
}

// 高性能流水线模式 - 双读取线程、实时处理、双输出线程
func optimizedPipelineMode(fq1, fq2, pattern, outdir string, percent int, numWorkers int, pigzPath string) error {
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

	// 创建高性能缓冲写入器
	bufferedW1 := NewBufferedFastqWriter(w1, outFile1, 8*1024*1024) // 8MB缓冲区
	bufferedW2 := NewBufferedFastqWriter(w2, outFile2, 8*1024*1024) // 8MB缓冲区
	defer bufferedW1.Close()
	defer bufferedW2.Close()

	// 创建scanner
	scanner1 := fastq.NewScanner(fastq.NewReader(reader1))
	scanner2 := fastq.NewScanner(fastq.NewReader(reader2))

	// 统计变量
	var totalReads, patternReads int64
	var keptPatternReads, totalOutputReads int
	var patternReadsInt int
	var mu sync.Mutex

	// 计算保留的pattern reads数量
	keepEveryN := 100 / percent

	// 创建写入缓存池 - 10M reads容量
	writeCache := make(chan *ProcessResult, 10000000) // 10M容量
	defer close(writeCache)

	// 创建读取队列 - 用于worker处理
	readQueue := make(chan ReadPair, 100000) // 10万容量
	defer close(readQueue)

	// 创建控制通道
	r2Queue := make(chan fastq.Sequence, 100000) // R2数据队列
	stopReading := make(chan bool)
	stopProcessing := make(chan bool)
	stopWriting := make(chan bool)
	r1BatchReady := make(chan []fastq.Sequence, 10) // R1批次就绪信号

	// 启动R1读取线程
	r1Done := make(chan bool)
	go func() {
		defer func() {
			r1Done <- true
		}()

		for scanner1.Next() {
			seq1 := scanner1.Seq()

			// 等待R2线程提供对应的seq2
			select {
			case r2Data := <-r2Queue:
				// 将read pair发送到处理队列
				select {
				case readQueue <- ReadPair{Seq1: seq1, Seq2: r2Data}:
					currentReads := atomic.AddInt64(&totalReads, 1)

					// 每读取2M reads输出进度
					if currentReads%2000000 == 0 {
						fmt.Printf("Read %d reads\n", currentReads)
					}
				case <-stopReading:
					return
				}
			case <-stopReading:
				return
			}
		}
	}()

	// 启动R2读取线程
	r2Done := make(chan bool)
	go func() {
		defer func() {
			r2Done <- true
		}()

		for scanner2.Next() {
			seq2 := scanner2.Seq()

			// 将R2数据发送到R1线程
			select {
			case r2Queue <- seq2:
			case <-stopReading:
				return
			}
		}
	}()

	// 启动worker线程池
	workerDone := make(chan bool)
	for i := 0; i < numWorkers; i++ {
		go func(workerID int) {
			defer func() {
				if workerID == 0 {
					workerDone <- true
				}
			}()

			for readPair := range readQueue {
				result := processReadPair(readPair.Seq1, readPair.Seq2, pattern, keepEveryN, &patternReadsInt, &mu)

				// 只有需要保留的reads才放入写入缓存池
				if result.Keep {
					select {
					case writeCache <- result:
					case <-stopProcessing:
						return
					}
				}
			}
		}(i)
	}

	// 启动R1输出线程
	r1WriteDone := make(chan bool)
	go func() {
		defer func() {
			r1WriteDone <- true
		}()

		batch := make([]fastq.Sequence, 0, 1000)

		for {
			select {
			case result := <-writeCache:
				if result == nil {
					// 结束信号
					if len(batch) > 0 {
						writeR1Batch(batch, bufferedW1)
					}
					return
				}

				batch = append(batch, result.Seq1)

				// 当达到1000条时，等待R2线程同步
				if len(batch) >= 1000 {
					select {
					case r1BatchReady <- batch:
						batch = make([]fastq.Sequence, 0, 1000)
					case <-stopWriting:
						return
					}
				}
			case <-stopWriting:
				return
			}
		}
	}()

	// 启动R2输出线程
	r2WriteDone := make(chan bool)
	go func() {
		defer func() {
			r2WriteDone <- true
		}()

		batch := make([]fastq.Sequence, 0, 1000)

		for {
			select {
			case result := <-writeCache:
				if result == nil {
					// 结束信号
					if len(batch) > 0 {
						writeR2Batch(batch, bufferedW2)
					}
					return
				}

				batch = append(batch, result.Seq2)

				// 当达到1000条时，等待R1线程同步
				if len(batch) >= 1000 {
					// 等待R1批次就绪
					select {
					case r1Batch := <-r1BatchReady:
						// 同步写入R1和R2
						writeR1Batch(r1Batch, bufferedW1)
						writeR2Batch(batch, bufferedW2)

						// 更新统计
						mu.Lock()
						totalOutputReads += len(batch)
						for _, r := range r1Batch {
							if containsPattern(string(r.Letters), pattern) {
								keptPatternReads++
							}
						}
						mu.Unlock()

						batch = make([]fastq.Sequence, 0, 1000)
					case <-stopWriting:
						return
					}
				}
			case <-stopWriting:
				return
			}
		}
	}()

	// 等待所有线程完成
	<-r1Done
	<-r2Done
	close(stopReading)
	close(stopProcessing)
	close(readQueue)

	// 发送结束信号给输出线程
	writeCache <- nil
	writeCache <- nil

	<-r1WriteDone
	<-r2WriteDone
	close(stopWriting)

	// 强制刷新缓冲区
	fmt.Println("Flushing all buffers...")
	if err := bufferedW1.Flush(); err != nil {
		fmt.Printf("Error flushing bufferedW1: %v\n", err)
	}
	if err := bufferedW2.Flush(); err != nil {
		fmt.Printf("Error flushing bufferedW2: %v\n", err)
	}

	// 检查错误
	if err := scanner1.Error(); err != nil {
		return fmt.Errorf("error reading fq1 (%s): %v", fq1, err)
	}
	if err := scanner2.Error(); err != nil {
		return fmt.Errorf("error reading fq2 (%s): %v", fq2, err)
	}

	// 计算最终比例
	finalRatio := 0.0
	if totalReads > 0 {
		finalRatio = float64(keptPatternReads) / float64(totalReads) * 100.0
	}

	// 输出统计信息
	statsFile := filepath.Join(outdir, basename2+".stats.txt")
	statsContent := fmt.Sprintf("Total reads: %d\nOriginal pattern reads: %d\nKept pattern reads: %d\nTotal output reads: %d\nFinal ratio: %.2f%%\n",
		totalReads, patternReads, keptPatternReads, totalOutputReads, finalRatio)

	if err := os.WriteFile(statsFile, []byte(statsContent), 0644); err != nil {
		return fmt.Errorf("failed to write stats file: %v", err)
	}

	// 输出最终结果
	fmt.Printf("\nOptimized pipeline mode completed. Results saved to: %s\n", statsFile)
	fmt.Printf("Filtered files saved to: %s and %s\n", outFile1Path, outFile2Path)
	fmt.Printf("Total reads: %d, Kept pattern reads: %d, Total output reads: %d, Final ratio: %.2f%%\n",
		totalReads, keptPatternReads, totalOutputReads, finalRatio)

	// 如果指定了pigz路径，进行压缩
	if pigzPath != "" {
		if err := compressWithPigz(outFile1Path, outFile2Path, pigzPath, numWorkers); err != nil {
			return fmt.Errorf("failed to compress output files: %v", err)
		}
	}

	return nil
}

// 内存映射的FastqReader
type MmapFastqReader struct {
	file      *os.File
	data      mmap.MMap
	scanner   *fastq.Scanner
	reader    io.Reader
	offset    int
	chunkSize int
	sequences []fastq.Sequence
	index     int
}

func NewMmapFastqReader(filePath string, chunkSize int) (*MmapFastqReader, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}

	var reader io.Reader = file
	var data mmap.MMap

	// 检查文件大小（用于调试）
	_, err = file.Stat()
	if err != nil {
		file.Close()
		return nil, err
	}

	// 如果文件是gzip格式，使用普通读取
	if strings.HasSuffix(filePath, ".gz") {
		gz, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, err
		}
		reader = gz
	} else {
		// 使用内存映射
		data, err = mmap.Map(file, mmap.RDONLY, 0)
		if err != nil {
			file.Close()
			return nil, err
		}
		reader = bytes.NewReader(data)
	}

	return &MmapFastqReader{
		file:      file,
		data:      data,
		reader:    reader,
		chunkSize: chunkSize,
		sequences: make([]fastq.Sequence, 0, 50000), // 预分配5万容量
		offset:    0,
		index:     0,
	}, nil
}

func (mfr *MmapFastqReader) Close() error {
	if mfr.data != nil {
		mfr.data.Unmap()
	}
	return mfr.file.Close()
}

func (mfr *MmapFastqReader) ReadChunk() ([]fastq.Sequence, error) {
	// 如果当前批次还有数据，直接返回
	if mfr.index < len(mfr.sequences) {
		start := mfr.index
		end := len(mfr.sequences)
		mfr.index = end
		return mfr.sequences[start:end], nil
	}

	// 如果使用内存映射
	if mfr.data != nil {
		end := mfr.offset + mfr.chunkSize
		if end > len(mfr.data) {
			end = len(mfr.data)
		}

		if mfr.offset >= len(mfr.data) {
			return nil, io.EOF
		}

		// 创建scanner处理这个块
		reader := bytes.NewReader(mfr.data[mfr.offset:end])
		scanner := fastq.NewScanner(fastq.NewReader(reader))

		// 清空之前的序列
		mfr.sequences = mfr.sequences[:0]
		mfr.index = 0

		// 解析FASTQ记录
		for scanner.Next() {
			mfr.sequences = append(mfr.sequences, scanner.Seq())
		}

		if err := scanner.Error(); err != nil {
			return nil, err
		}

		mfr.offset = end
		return mfr.sequences, nil
	}

	// 如果使用普通读取（gzip文件）
	buffer := make([]byte, mfr.chunkSize)
	n, err := mfr.reader.Read(buffer)
	if err != nil && err != io.EOF {
		return nil, err
	}
	if n == 0 {
		return nil, io.EOF
	}

	// 创建scanner处理这个块
	reader := bytes.NewReader(buffer[:n])
	scanner := fastq.NewScanner(fastq.NewReader(reader))

	// 清空之前的序列
	mfr.sequences = mfr.sequences[:0]
	mfr.index = 0

	// 解析FASTQ记录
	for scanner.Next() {
		mfr.sequences = append(mfr.sequences, scanner.Seq())
	}

	if err := scanner.Error(); err != nil {
		return nil, err
	}

	return mfr.sequences, nil
}

// 内存映射高性能流水线模式
func mmapPipelineMode(fq1, fq2, pattern, outdir string, percent int, numWorkers int, pigzPath string) error {
	// 创建输出目录
	if err := os.MkdirAll(outdir, 0755); err != nil {
		return fmt.Errorf("failed to create output directory: %v", err)
	}

	// 创建内存映射读取器，5MB缓冲区
	chunkSize := 5 * 1024 * 1024 // 5MB块
	reader1, err := NewMmapFastqReader(fq1, chunkSize)
	if err != nil {
		return err
	}
	defer reader1.Close()

	reader2, err := NewMmapFastqReader(fq2, chunkSize)
	if err != nil {
		return err
	}
	defer reader2.Close()

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
	bufferedW1 := NewBufferedFastqWriter(w1, outFile1, 16*1024*1024) // 16MB缓冲区
	bufferedW2 := NewBufferedFastqWriter(w2, outFile2, 16*1024*1024) // 16MB缓冲区
	defer bufferedW1.Close()
	defer bufferedW2.Close()

	// 统计变量
	var totalReads, patternReads int64
	var keptPatternReads, totalOutputReads int
	var patternReadsInt int
	var mu sync.Mutex

	// 计算保留的pattern reads数量
	keepEveryN := 100 / percent

	// 创建写入缓存池，增加容量
	writeCache := make(chan *ProcessResult, 20000000) // 20M容量
	defer close(writeCache)

	// 创建读取队列，增加容量
	readQueue := make(chan ReadPair, 500000) // 50万容量
	defer close(readQueue)

	// 创建控制通道
	stopReading := make(chan bool)
	stopProcessing := make(chan bool)
	stopWriting := make(chan bool)
	r1BatchReady := make(chan []fastq.Sequence, 20) // 增加缓冲区

	// 启动多个读取线程（4个线程）
	numReadThreads := 4
	readDone := make(chan bool, numReadThreads)

	for threadID := 0; threadID < numReadThreads; threadID++ {
		go func(threadID int) {
			defer func() {
				readDone <- true
			}()

			for {
				// 读取R1块
				r1Sequences, err := reader1.ReadChunk()
				if err == io.EOF {
					break
				}
				if err != nil {
					log.Printf("Error reading R1 chunk in thread %d: %v", threadID, err)
					break
				}

				// 读取R2块
				r2Sequences, err := reader2.ReadChunk()
				if err == io.EOF {
					break
				}
				if err != nil {
					log.Printf("Error reading R2 chunk in thread %d: %v", threadID, err)
					break
				}

				// 处理这个块中的reads
				minLen := len(r1Sequences)
				if len(r2Sequences) < minLen {
					minLen = len(r2Sequences)
				}

				for i := 0; i < minLen; i++ {
					select {
					case readQueue <- ReadPair{Seq1: r1Sequences[i], Seq2: r2Sequences[i]}:
						currentReads := atomic.AddInt64(&totalReads, 1)

						// 每读取2M reads输出进度
						if currentReads%2000000 == 0 {
							fmt.Printf("Read %d reads (mmap mode, thread %d)\n", currentReads, threadID)
						}
					case <-stopReading:
						return
					}
				}
			}
		}(threadID)
	}

	// 启动worker线程池
	workerDone := make(chan bool)
	for i := 0; i < numWorkers; i++ {
		go func(workerID int) {
			defer func() {
				if workerID == 0 {
					workerDone <- true
				}
			}()

			for readPair := range readQueue {
				result := processReadPair(readPair.Seq1, readPair.Seq2, pattern, keepEveryN, &patternReadsInt, &mu)

				if result.Keep {
					select {
					case writeCache <- result:
					case <-stopProcessing:
						return
					}
				}
			}
		}(i)
	}

	// 启动输出线程
	r1WriteDone := make(chan bool)
	go func() {
		defer func() {
			r1WriteDone <- true
		}()

		batch := make([]fastq.Sequence, 0, 1000)

		for {
			select {
			case result := <-writeCache:
				if result == nil {
					if len(batch) > 0 {
						writeR1Batch(batch, bufferedW1)
					}
					return
				}

				batch = append(batch, result.Seq1)

				if len(batch) >= 1000 {
					select {
					case r1BatchReady <- batch:
						batch = make([]fastq.Sequence, 0, 1000)
					case <-stopWriting:
						return
					}
				}
			case <-stopWriting:
				return
			}
		}
	}()

	r2WriteDone := make(chan bool)
	go func() {
		defer func() {
			r2WriteDone <- true
		}()

		batch := make([]fastq.Sequence, 0, 1000)

		for {
			select {
			case result := <-writeCache:
				if result == nil {
					if len(batch) > 0 {
						writeR2Batch(batch, bufferedW2)
					}
					return
				}

				batch = append(batch, result.Seq2)

				if len(batch) >= 1000 {
					select {
					case r1Batch := <-r1BatchReady:
						writeR1Batch(r1Batch, bufferedW1)
						writeR2Batch(batch, bufferedW2)

						mu.Lock()
						totalOutputReads += len(batch)
						for _, r := range r1Batch {
							if containsPattern(string(r.Letters), pattern) {
								keptPatternReads++
							}
						}
						mu.Unlock()

						batch = make([]fastq.Sequence, 0, 1000)
					case <-stopWriting:
						return
					}
				}
			case <-stopWriting:
				return
			}
		}
	}()

	// 等待所有读取线程完成
	for i := 0; i < numReadThreads; i++ {
		<-readDone
	}
	close(stopReading)
	close(stopProcessing)
	close(readQueue)

	writeCache <- nil
	writeCache <- nil

	<-r1WriteDone
	<-r2WriteDone
	close(stopWriting)

	// 强制刷新缓冲区
	fmt.Println("Flushing all buffers...")
	if err := bufferedW1.Flush(); err != nil {
		fmt.Printf("Error flushing bufferedW1: %v\n", err)
	}
	if err := bufferedW2.Flush(); err != nil {
		fmt.Printf("Error flushing bufferedW2: %v\n", err)
	}

	// 计算最终比例
	finalRatio := 0.0
	if totalReads > 0 {
		finalRatio = float64(keptPatternReads) / float64(totalReads) * 100.0
	}

	// 输出统计信息
	statsFile := filepath.Join(outdir, basename2+".stats.txt")
	statsContent := fmt.Sprintf("Total reads: %d\nOriginal pattern reads: %d\nKept pattern reads: %d\nTotal output reads: %d\nFinal ratio: %.2f%%\n",
		totalReads, patternReads, keptPatternReads, totalOutputReads, finalRatio)

	if err := os.WriteFile(statsFile, []byte(statsContent), 0644); err != nil {
		return fmt.Errorf("failed to write stats file: %v", err)
	}

	// 输出最终结果
	fmt.Printf("\nMemory-mapped pipeline mode completed. Results saved to: %s\n", statsFile)
	fmt.Printf("Filtered files saved to: %s and %s\n", outFile1Path, outFile2Path)
	fmt.Printf("Total reads: %d, Kept pattern reads: %d, Total output reads: %d, Final ratio: %.2f%%\n",
		totalReads, keptPatternReads, totalOutputReads, finalRatio)

	// 如果指定了pigz路径，进行压缩
	if pigzPath != "" {
		if err := compressWithPigz(outFile1Path, outFile2Path, pigzPath, numWorkers); err != nil {
			return fmt.Errorf("failed to compress output files: %v", err)
		}
	}

	return nil
}

// 写入R1批次
func writeR1Batch(batch []fastq.Sequence, writer *BufferedFastqWriter) {
	for _, seq := range batch {
		if err := writer.Write(seq); err != nil {
			log.Printf("Error writing R1: %v", err)
		}
	}
}

// 写入R2批次
func writeR2Batch(batch []fastq.Sequence, writer *BufferedFastqWriter) {
	for _, seq := range batch {
		if err := writer.Write(seq); err != nil {
			log.Printf("Error writing R2: %v", err)
		}
	}
}

func main() {
	// 定义命令行参数
	var fq1, fq2, pattern, outdir, pigzPath string
	var percent, numWorkers int
	var showVersion bool
	var useMmap bool

	flag.StringVar(&fq1, "fq1", "", "Input fastq file 1 (supports .gz format)")
	flag.StringVar(&fq2, "fq2", "", "Input fastq file 2 (supports .gz format)")
	flag.StringVar(&pattern, "pattern", "AGCAGTGGTATCAACGCAGAGTACA", "Pattern to search for")
	flag.StringVar(&outdir, "outdir", "", "Output directory")
	flag.IntVar(&percent, "percent", 5, "Percentage of pattern reads to keep (0-100)")
	flag.IntVar(&numWorkers, "workers", 4, "Number of worker threads")
	flag.StringVar(&pigzPath, "pigz", "", "Path to pigz executable for compression")
	flag.BoolVar(&showVersion, "version", false, "Show version information")
	flag.BoolVar(&useMmap, "mmap", false, "Use memory-mapped file reading for better performance")

	// 解析命令行参数
	flag.Parse()

	// 显示版本信息
	if showVersion {
		fmt.Printf("patternqc v%s (Build: %s)\n", Version, BuildTime)
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

	// 根据参数选择执行模式
	var err error
	if useMmap {
		fmt.Println("Using memory-mapped pipeline mode...")
		err = mmapPipelineMode(fq1, fq2, pattern, outdir, percent, numWorkers, pigzPath)
	} else {
		fmt.Println("Using optimized pipeline mode...")
		err = optimizedPipelineMode(fq1, fq2, pattern, outdir, percent, numWorkers, pigzPath)
	}
	if err != nil {
		log.Fatalf("Error: %v", err)
	}
}
