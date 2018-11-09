package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"os"
	"strings"

	"github.com/jessevdk/go-flags"
)

func computeDistance(filename string, threshold int, nColumns int, separator string) {

	//Prints a sparse matrix containing all pair-wise similarities between the input lines to stdout

	myFile, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer myFile.Close()

	record := csv.NewReader(myFile)
	if separator == "" {
		record.Comma = '\t'
	} else {
		record.Comma = []rune(separator)[0]
	}
	lines, err := record.ReadAll()

	if err != nil {
		panic(err)
	}

	for i, lineX := range lines {
		for j, lineY := range lines[i+1:] {
			matchingCount := matchingDistance(lineX, lineY, nColumns)
			if matchingCount >= threshold {
				fmt.Println(i+1, j+1, matchingCount)
			}
		}
	}
}

func findCores(filename string) {

	myFile, err := os.Open(filename)
	defer myFile.Close()

	if err != nil {
		panic(err)
	}

	encountered := map[string][]int{}

	scanner := bufio.NewScanner(myFile)

	index := 1
	for scanner.Scan() {
		currentLine := scanner.Text()
		_, wasEncountered := encountered[currentLine]
		if wasEncountered {
			encountered[currentLine] = append(encountered[currentLine], index)
		} else {
			newIndex := []int{index}
			encountered[currentLine] = newIndex
		}
		index++
	}

	for _, indices := range encountered {
		lineToWrite := strings.Trim(fmt.Sprint(indices), "[]")

		fmt.Println(len(indices), lineToWrite)
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}

}

func matchingDistance(u []string, v []string, nColumns int) int {

	//Compute the Hamming similarity between two non-binary vectors

	n := len(u)
	m := len(v)

	if n != m {
		panic("Row lengths don't match - uneven number of fields in input file")
	}

	nColumnsToConsider := n

	if nColumns <= n && nColumns > 0 {
		nColumnsToConsider = nColumns
	}

	matchingCount := 0
	for i := 0; i < nColumnsToConsider; i++ {
		if u[i] == v[i] {
			matchingCount++
		}
	}

	return matchingCount
}

func main() {

	//A very rudimentary parser for the time being

	type Options struct {
		Input string `short:"i" long:"input" description:"Input community file in csv format" required:"true"`

		Threshold int `short:"t" long:"threshold" description:"Minimum matches threshold" default:"1"`

		Separator string `short:"s" long:"separator" description:"Field separator for csv input" default:""`

		nColumns int `short:"c" long:"columns" description:"Number of columns to consider" default:"0"`

		Mode string `short:"m" long:"mode" description:"Mode" required:"true"`
	}

	var options Options

	_, err := flags.Parse(&options)

	if err != nil {
		fmt.Println("Error when parsing arguments")
		panic(err)
	}

	if options.Mode == "cores" {
		findCores(options.Input)
	} else if options.Mode == "distances" {
		threshold := options.Threshold
		nColumns := options.nColumns
		separator := options.Separator
		computeDistance(options.Input, threshold, nColumns, separator)
	} else {
		panic("Wrong mode. Available modes are: distances, cores")
	}
}
