/*
count k-mer frequencies in a set of reads, using multiple goroutines.
*/

package main

import (
   "github.com/vtphan/kmers"
   "os"
   "bufio"
   // "fmt"
   "runtime"
   "sync"
   "math"
   "encoding/csv"
   "strconv"
   "io/ioutil"
)


func CountFreq(readFile string, freq map[int]int, K int) {

   // Get all reads into a channel
   reads := make(chan string)
   go func() {
      f, err := os.Open(readFile)
      if err != nil {
         panic("Error opening " + readFile)
      }
      defer f.Close()
      scanner := bufio.NewScanner(f)
      i := 0.0
      for scanner.Scan() {
        if math.Mod(i, 2.0) == 1 {
         reads <- (scanner.Text())
        }
        i++
      }
      close(reads)
   }()

   numCores := runtime.NumCPU()
   runtime.GOMAXPROCS(numCores)
   var wg sync.WaitGroup

   // Here is a map of kmers that need to be counted.  Only these kmers are counted.
/*   freq := make(map[int]int)
   freq[158] = 0
   freq[180] = 0
   freq[39] = 0
   freq[59] = 0*/

   // Start a new counter that counts only kmers in freq.
   c := kmers.NewKmerCounter(K, freq)

   // Count kmers in different cores simultaneously.
   for i:=0; i<numCores; i++ {
      wg.Add(1)
      go func() {
         defer wg.Done()
         for read := range(reads){
            // Count2 counts on both strands. Count1 counts only on the main strand.
            c.Count1([]byte(read))
         }
      }()
   }

   // Finish counting
   wg.Wait()

   // Print out the result
/*   for kmer := range(freq) {
      fmt.Println(kmer,kmers.NumToKmer(kmer,K),freq[kmer])
   }*/
}
func error_check(err error) {
   if err != nil {
      panic(err)
   }
}
func getGSM(GSMfile string, freq map[int]int) {
   csvfile, err := os.Open(GSMfile)
   error_check(err)
   defer csvfile.Close()
   reader := csv.NewReader(csvfile)
   _, err = reader.Read()

   for {
      line, err := reader.Read()
      if err != nil {
          break
      } else {
          kmer, _ := strconv.Atoi(line[0])
          freq[kmer] = 0
      }
   }
}
func assignTocsv(GSMfile string, freq map[int]int) {
   resultfile, err := os.Create("b/" + GSMfile)
   error_check(err)
   rw := csv.NewWriter(resultfile)
   head := make([]string, 2)
   head[0] = "kmer"
   head[1] = "b"
   returnError := rw.Write(head)
   error_check(returnError)
   rw.Flush()

   csvfile, err := os.Open("genomeGSM/" + GSMfile)
   error_check(err)
   defer csvfile.Close()
   reader := csv.NewReader(csvfile)

   for {
      line, err := reader.Read()
      if err != nil {
          break
      } else {
         head[0] = line[0]
         kmer, _ := strconv.Atoi(line[0])
         // fmt.Println(kmers.NumToKmer(kmer,14))
         head[1] = strconv.Itoa(freq[kmer])
         returnError := rw.Write(head)
         error_check(returnError)
      }
   }
   rw.Flush()
}
func main() {
   if len(os.Args) != 5 {
      panic("Must provide readfile and GSM, GSM folder from genomes, kmer length!")
   }
   _ = os.Mkdir("b", 0777)
   files, _ := ioutil.ReadDir(os.Args[3])
   K, _ := strconv.Atoi(os.Args[4])
   freq := make(map[int]int)
   getGSM(os.Args[2], freq)
   CountFreq(os.Args[1], freq, K)
   for _, file := range files{
      assignTocsv(file.Name(), freq)
   }
}
