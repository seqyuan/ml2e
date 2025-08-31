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

	"github.com/seqyuan/annogene/io/fastq"
)

const (
	Version   = "3.0.3"
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

// 简化的流水线模式 - 10M reads pair一个单位，处理完立即写入
func simplePipelineMode(fq1, fq2, pattern, outdir string, percent int, numWorkers int, pigzPath string) error {
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

	// 配置参数
	batchSize := 10000000 // 10M reads per batch

	// 统计变量
	var totalReads, patternReads, keptPatternReads, totalOutputReads int64
	var patternReadsInt int // 用于processReadPair函数
	var mu sync.Mutex

	// 计算保留的pattern reads数量
	keepEveryN := 100 / percent

	// 主处理循环
	batchID := 0
	for {
		// 读取一个批次
		readPairs := make([]ReadPair, 0, batchSize)

		// 读取10M reads或直到文件结束
		for i := 0; i < batchSize; i++ {
			if !scanner1.Next() || !scanner2.Next() {
				break
			}

			seq1 := scanner1.Seq()
			seq2 := scanner2.Seq()

			readPairs = append(readPairs, ReadPair{
				Seq1: seq1,
				Seq2: seq2,
			})

			atomic.AddInt64(&totalReads, 1)
		}

		// 如果没有读取到数据，退出
		if len(readPairs) == 0 {
			break
		}

		fmt.Printf("Processing batch %d: %d reads, Total: %d\n", batchID, len(readPairs), totalReads)

		// 处理批次 - 使用多个worker并行处理
		results := make([]*ProcessResult, len(readPairs))
		var wg sync.WaitGroup

		// 将批次分成多个子批次给worker处理
		subBatchSize := len(readPairs) / numWorkers
		if subBatchSize == 0 {
			subBatchSize = 1
		}

		for i := 0; i < numWorkers; i++ {
			wg.Add(1)
			go func(workerID int) {
				defer wg.Done()

				start := workerID * subBatchSize
				end := start + subBatchSize
				if workerID == numWorkers-1 {
					end = len(readPairs) // 最后一个worker处理剩余的所有
				}

				for j := start; j < end; j++ {
					result := processReadPair(readPairs[j].Seq1, readPairs[j].Seq2, pattern, keepEveryN, &patternReadsInt, &mu)
					results[j] = result
				}
			}(i)
		}

		// 等待所有worker完成
		wg.Wait()

		// 立即写入结果
		r1Data := make([]fastq.Sequence, 0, len(results))
		r2Data := make([]fastq.Sequence, 0, len(results))

		for _, result := range results {
			if result.Keep {
				r1Data = append(r1Data, result.Seq1)
				r2Data = append(r2Data, result.Seq2)
			}
		}

		// 并行写入R1和R2
		var writeWg sync.WaitGroup
		writeWg.Add(2)

		// 写入R1
		go func() {
			defer writeWg.Done()
			for _, seq := range r1Data {
				if err := bufferedW1.Write(seq); err != nil {
					log.Printf("Error writing R1: %v", err)
				}
			}
		}()

		// 写入R2
		go func() {
			defer writeWg.Done()
			for _, seq := range r2Data {
				if err := bufferedW2.Write(seq); err != nil {
					log.Printf("Error writing R2: %v", err)
				}
			}
		}()

		// 等待写入完成
		writeWg.Wait()

		// 更新统计
		mu.Lock()
		totalOutputReads += int64(len(r1Data))
		for _, result := range results {
			if result.Keep && result.HasPattern {
				keptPatternReads++
			}
		}
		mu.Unlock()

		fmt.Printf("Batch %d completed: %d reads written, Total output: %d\n", batchID, len(r1Data), totalOutputReads)

		batchID++

		// 如果读取的数据少于批次大小，说明文件结束了
		if len(readPairs) < batchSize {
			break
		}
	}

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

	// 输出统计信息
	statsFile := filepath.Join(outdir, basename2+".stats.txt")
	statsContent := fmt.Sprintf("Total reads: %d\nOriginal pattern reads: %d\nKept pattern reads: %d\nTotal output reads: %d\nFinal ratio: %.2f%%\n",
		totalReads, patternReads, keptPatternReads, totalOutputReads, finalRatio)

	if err := os.WriteFile(statsFile, []byte(statsContent), 0644); err != nil {
		return fmt.Errorf("failed to write stats file: %v", err)
	}

	// 输出最终结果
	fmt.Printf("\nFilter mode completed. Results saved to: %s\n", statsFile)
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
	err := simplePipelineMode(fq1, fq2, pattern, outdir, percent, numWorkers, pigzPath)
	if err != nil {
		log.Fatalf("Error: %v", err)
	}
}
