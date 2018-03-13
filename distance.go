package main

import (
	"bufio"
	"encoding/csv"
	"github.com/jessevdk/go-flags"
	"fmt"
	"os"
	"strings"
)

func compute_distance(filename string, threshold int, n_columns int, separator string) {

	//Prints a sparse matrix containing all pair-wise similarities between the input lines to stdout
	
	my_file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer my_file.Close()

	record := csv.NewReader(my_file)
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
		for j, lineY := range lines {
			matching_count := matching_distance(lineX, lineY, n_columns)
			if matching_count >= threshold {
				fmt.Println(i+1, j+1, matching_count)
			}
		}
	}
}

func find_cores(filename string) {

	my_file, err := os.Open(filename)
	defer my_file.Close()

	if err != nil {
		panic(err)
	}

	encountered := map[string][]int{}

	scanner := bufio.NewScanner(my_file)

	index := 1
	for scanner.Scan() {
		current_line := scanner.Text()
		_, was_encountered := encountered[current_line]
		if was_encountered {
			encountered[current_line] = append(encountered[current_line], index)
		} else {
			new_index := []int{index}
			encountered[current_line] = new_index
		}
		index += 1
	}

	for _, indices := range encountered {
		line_to_write := strings.Trim(fmt.Sprint(indices), "[]")

		fmt.Println(len(indices), line_to_write)
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}

}

func matching_distance(u []string, v []string, n_columns int) int {

	//Compute the Hamming similarity between two non-binary vectors

	n := len(u)
	m := len(v)

	if n != m {
		panic("Row lengths don't match - uneven number of fields in input file")
	}

	n_columns_to_consider := n

	if n_columns <= n && n_columns > 0 {
		n_columns_to_consider = n_columns
	}

	matching_count := 0
	for i := 1; i < n_columns_to_consider; i++ {
		if u[i] == v[i] {
			matching_count += 1
		}
	}

	return matching_count
}

func main() {

	//A very rudimentary parser for the time being

	type Options struct {

	Input string `short:"i" long:"input" description:"Input community file in csv format" required:"true"`

	Threshold int `short:"t" long:"threshold" description:"Minimum matches threshold" default:"1"`

	Separator string `short:"s" long:"separator" description:"Field separator for csv input" default:""`

	N_columns int `short:"c" long:"columns" description:"Number of columns to consider" default:"0"`

	Mode string `short:"m" long:"mode" description:"Mode" required:"true"`
	}

	var options Options

	_, err := flags.Parse(&options)

	if err != nil {
		fmt.Println("Error when parsing arguments")
		panic(err)
	}

	if options.Mode == "cores" {
		find_cores(options.Input)
	} else if options.Mode == "distances" {
		threshold := options.Threshold
		n_columns := options.N_columns
		separator := options.Separator
		compute_distance(options.Input, threshold, n_columns, separator)
	} else {
		panic("Wrong mode. Available modes are: distances, cores")
	}
}
