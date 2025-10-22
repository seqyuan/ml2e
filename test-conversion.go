package main

import (
	"fmt"
	"strings"
)

// 修改readID中的@ML15为@E251
func modifyReadID(readID string) string {
	if strings.HasPrefix(readID, "@ML15") {
		return "@E251" + readID[5:] // 将@ML15替换为@E251
	}
	return readID
}

func main() {
	// 测试用例
	testCases := []string{
		"@ML1512345:1:1101:1000:2211/1",
		"@ML1512345:1:1101:1000:2211/2",
		"@AB12345:1:1101:1000:2211/1",
		"@ML12345:1:1101:1000:2211/1",
		"@ML15",
		"@ML15ABC",
	}

	fmt.Println("🧪 测试 readID 转换功能：")
	fmt.Println("==================================================")

	for _, testCase := range testCases {
		result := modifyReadID(testCase)
		fmt.Printf("输入: %-30s → 输出: %s\n", testCase, result)
	}

	fmt.Println("\n✅ 测试完成！")
}
