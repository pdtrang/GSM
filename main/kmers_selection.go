/*
   Count k-mers from FASTA files efficiently in one pass, using O(1) space.
*/
package main

import (
   "os"
   "io"
   "io/ioutil"
   "path"
   "strconv"
   "math"
   "sync"
   "runtime"
   "fmt"
   "strings"
   "sort"
   "encoding/csv"
)

var BUFF_SIZE = 1024 * 1024   // default buffer size is 1MB
const MaxInt = int(^uint(0)>>1)
var MaxK = int(math.Log2(float64(MaxInt)) / 2.0)
var Code = map[byte]int { 'A':0, 'C':1, 'G':2, 'T':3 }
var CodeRC = map[byte]int { 'A':3, 'C':2, 'G':1, 'T': 0 }
var Neighbor = [][]int{}

type ByArray struct {
   weights []int
   keys []int32
}
func (s ByArray) Len() int {
   return len(s.keys)
}
func (s ByArray) Swap(i ,j int) {
   s.keys[i], s.keys[j] = s.keys[j], s.keys[i]
}
func (s ByArray) Less(i ,j int) bool {
   return (s.weights[s.keys[i]] < s.weights[s.keys[j]])
}

func SetBufferSize(size int) {
   BUFF_SIZE = size
}

func error_check(err error) {
   if err != nil {
      panic(err)
   }
}

func AssignKmersToGenomes(occ []int16, K int, N int) [][]int {
   genomes := make([][]int, N)
   for i:=0; i<len(occ); i++ {
      if occ[i] > 0 {
         // fmt.Printf("%10d %3d %s\n", i, occ[i], NumToKmer(i, K))
         genomes[occ[i]-1] = append(genomes[occ[i]-1],i)
      }
   }
   return genomes
}

func CountUniqueKmers(dirname string, K int) ([]int16, []int, []int, int) {
   var wg sync.WaitGroup
   occ := make([]int16, 1<<uint(2*K))
   pos := make([]int, 1<<uint(2*K))
   freq := make([]int, 1<<uint(2*K))
   lock := make([]sync.Mutex, 1<<uint(2*K))
   files, err := ioutil.ReadDir(dirname)
   error_check(err)

   wg.Add(len(files))
   namesChan := make(chan string, len(files))
   for idx, f := range(files) {
      namesChan <- (strconv.Itoa(idx) + "," + f.Name())
   }
   close(namesChan)

   workers := runtime.NumCPU()
   if len(files) < workers {
      workers = len(files)
   }

   for i := 0; i < workers; i++ {
      go func(){
         for name := range namesChan {
            elements := strings.Split(name, ",")
            idx, _ := strconv.Atoi(elements[0])
            filename := elements[1]
            // fmt.Println("processing", filename)
            file, err := os.Open(path.Join(dirname, filename))
            error_check(err)
            countFASTA(file, idx+1, K, occ, pos, freq, lock)
            file.Close()
            wg.Done()
         }
      }()
   }

   wg.Wait()
   return occ, pos, freq, len(files)
}

func countFASTA(file *os.File, gid int, K int, occ []int16, pos []int, freq []int, lock []sync.Mutex) {
   buffer := make([]byte, BUFF_SIZE)
   previous_kmer := make([]byte, K)
   header := true
   n, j, value, value_rc := 0, 0, 0, 0
   var pk int
   var err error
   // base * 4^K
   fourK := map[byte]int {'A': 0, 'C': 1<<(uint(2*K)), 'G': 2<<(uint(2*K)), 'T': 3<<(uint(2*K)) }
   fourKRC := map[byte]int { 'A':fourK['T'], 'C':fourK['G'], 'G':fourK['C'], 'T':fourK['A'] }
   k := 0

   for {
      n, err = file.Read(buffer)
      if err != nil {
         if err == io.EOF {
            break;
         }
         panic(err)
      }
      // fmt.Println("buffer:", string(buffer[0:n]))
      for i:=0; i<n; i++ {
         if header {
            if buffer[i] == '\n' {
               header = false
               // fmt.Println("new contig starting at", string(buffer[i+1:n]))
            }
         } else {
            if buffer[i] == '>' {
               header = true
               j, value, value_rc = 0, 0, 0
            } else if buffer[i] == '\n' {
               ;
            } else if buffer[i] == 'N' {
               j, value, value_rc = 0, 0, 0
            } else {
               if j < K {
                  previous_kmer[pk] = buffer[i]
                  value = value<<2 + Code[buffer[i]]
                  value_rc = value_rc + CodeRC[buffer[i]]<<uint(j<<1)
                  j++
               } else {
                  value = value<<2 + Code[buffer[i]] - fourK[previous_kmer[pk]]
                  value_rc = (value_rc + fourKRC[buffer[i]] - CodeRC[previous_kmer[pk]]) >> 2
                  previous_kmer[pk] = buffer[i]
               }
               pk = (pk+1)%len(previous_kmer)

               if j==K {
                  lock[value].Lock()
                  freq[value]++
                  if occ[value] == 0 {
                     occ[value] = int16(gid)
                     pos[value] = i + BUFF_SIZE * k
                  } else if occ[value] > 0 && occ[value] != int16(gid) {
                     occ[value] = int16(-1)      // value is no longer a GSM
                  }
                  lock[value].Unlock()

                  lock[value_rc].Lock()
                  freq[value_rc]++
                  if occ[value_rc] == 0 {
                     occ[value_rc] = int16(gid)
                     pos[value_rc] = i + BUFF_SIZE * k
                  } else if occ[value_rc] > 0 && occ[value_rc] != int16(gid) {
                     occ[value_rc] = -1   // value_rc is no longer a GSM
                  }
                  lock[value_rc].Unlock()
               }
            }
         }
      }
      k++
   }
}

func NumToKmer(x int, K int) string {
   y := make([]byte, K)
   for i:=K-1; i>=0; i-- {
      base := x%4
      switch base {
         case 0: y[i] = 'A'
         case 1: y[i] = 'C'
         case 2: y[i] = 'G'
         case 3: y[i] = 'T'
      }
      x = (x-base)>>2
   }
   return string(y)
}

func InitNeighbors(K int) {
   Neighbor = make([][]int, K)
   for i:=0; i<len(Neighbor); i++ {
      Neighbor[i] = make([]int, 4)
      for j:=1; j<4; j++ {
         Neighbor[i][j] = j << uint(2*(K-i-1))
      }
   }
   for i:=0; i<len(Neighbor); i++ {
      fmt.Println(Neighbor[i])
   }
}

func GetNeighbors(x, K int) []int{
   digits := make([]int, K)
   n1 := make([]int, 3*K)
   v := x
   var pow uint

   for i:=K-1; i>=0; i-- {
      digits[i] = x%4
      x = (x-digits[i])>>2
   }
   for i, j:=0, 0; i<K; i++ {
      pow = uint(2*(K-1-i))
      switch digits[i] {
      case 0:
         n1[j] = v + 1 << pow
         n1[j+1] = v + 2 << pow
         n1[j+2] = v + 3 << pow
      case 1:
         n1[j] = v - 1 << pow
         n1[j+1] = v + 1 << pow
         n1[j+2] = v + 2 << pow
      case 2:
         n1[j] = v - 2 << pow
         n1[j+1] = v - 1 << pow
         n1[j+2] = v + 1 << pow
      case 3:
         n1[j] = v - 3 << pow
         n1[j+1] = v - 2 << pow
         n1[j+2] = v - 1 << pow
      }
      j += 3
   }

   return n1
}

func main() {
   runtime.GOMAXPROCS(runtime.NumCPU())
   
   if len(os.Args) != 4 {
      println("\tgo run kmers_selection.go genomes_dir K goodnum")
      println("\t\tgenomes_dir : directory containing genomes in FASTA format.")
      println("\t\tK : kmer length.  Should not be more than 16.")
      println("\t\tgoodnum : number of good kmers you need.")
      os.Exit(1)
   }
   K, _ := strconv.Atoi(os.Args[2])
   goodnum, _ := strconv.Atoi(os.Args[3])
   
   if K >= MaxK {
      println(K,"is too big.")
      os.Exit(1)
   }
   
   // InitNeighbors(K)
   gsmfile, err := os.Create("gsm.csv")
   if err != nil {
      fmt.Printf("%v\n", err)
      os.Exit(1)
   }
   gsm := csv.NewWriter(gsmfile)
   defer gsmfile.Close()
   files, _ := ioutil.ReadDir(os.Args[1])
   head := make([]string, len(files)+1)
   head0 := make([]string, 2)
   headpos := make([]string, 3)
   head[0] = "kmer"
   for index, fi := range files {
      head[index+1] = fi.Name()
   }
   returnError := gsm.Write(head)
   if returnError != nil {
      fmt.Println(returnError)
   }
   gsm.Flush()

   posfile, err := os.Create("pos.csv")
   if err != nil {
      fmt.Printf("%v\n", err)
      os.Exit(1)
   }
   position := csv.NewWriter(posfile)
   defer posfile.Close()

   occ, pos, freq, N := CountUniqueKmers(os.Args[1], K)
   genomes := AssignKmersToGenomes(occ, K, N)
   var wg sync.WaitGroup
   wg.Add(len(genomes))
   namesChan := make(chan int, len(genomes))
   for i := 0; i < len(genomes); i++ {
      namesChan <- i
   }
   close(namesChan)
   workers := 1
   if len(genomes) < workers {
      workers = len(genomes)
   }
   _ = os.Mkdir("genomeGSM", 0777)

   for i := 0; i < workers; i++ {
      go func(){
         for i := range namesChan {
            gsm0file, _ := os.Create("genomeGSM/" + files[i].Name() + ".csv")
            gsm0 := csv.NewWriter(gsm0file)
            // fmt.Println(i, len(genomes[i]))
            weights := make([]int, len(genomes[i]))
            keys := make([]int32, len(genomes[i]))
            for j:=0; j<len(genomes[i]); j++ {
               keys[j] = int32(j)
               // fmt.Println("\tGSM", NumToKmer(genomes[i][j], K), genomes[i][j])
               neighbors := GetNeighbors(genomes[i][j], K)
               for m := 0; m < len(neighbors); m++ {
                  weights[j] = weights[j] + freq[neighbors[m]]
               }
            }
            sort.Sort(ByArray{weights, keys})
            // fmt.Println("weights", weights)
            // fmt.Println(keys)
            if len(keys) != 0 {
               for n := 0; n < int(math.Min(float64(goodnum), float64(len(genomes[i])))); n++ {
                  // fmt.Println(keys[n], genomes[i][keys[n]])
                  // fmt.Println(genomes[i][keys[n]], weights[keys[n]])
                  head[0] = strconv.Itoa(genomes[i][keys[n]])
                  head0[0] = strconv.Itoa(genomes[i][keys[n]])
                  for j := 1; j <= len(files); j++ {
                     head[j] = strconv.Itoa(0)
                  }
                  head[i+1] = strconv.Itoa(freq[genomes[i][keys[n]]])
                  head0[1] = strconv.Itoa(freq[genomes[i][keys[n]]])
                  headpos[0] = strconv.Itoa(genomes[i][keys[n]])
                  headpos[1] = strconv.Itoa(int(occ[genomes[i][keys[n]]]))
                  headpos[2] = strconv.Itoa(pos[genomes[i][keys[n]]])
                  returnError := gsm.Write(head)
                  _ = gsm0.Write(head0)
                  _ = position.Write(headpos)
                  if returnError != nil {
                     fmt.Println(returnError)
                  }
                  gsm.Flush()
                  position.Flush()
                  // fmt.Println(i, NumToKmer(genomes[i][keys[n]], K), genomes[i][keys[n]], weights[keys[n]])
               }
               gsm0.Flush()
               gsm0file.Close()
               // fmt.Println("MaxWeight", weights[keys[int(math.Min(float64(goodnum), float64(len(genomes[i])))) - 1]])
            }
            wg.Done()
         }
      }()
   }
   wg.Wait()
   
}
