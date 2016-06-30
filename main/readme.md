##Usage

### Selecting GSMs

```
go run kmers_selection.go refgenomes_dir K numberOfGSMs
```
**Input**
- refgenomes_dir: contains Reference Genomes.
- K: length of marker.
- numberOfGSMs: number of markers that the user want to select. Markers are selected by lowest global and highest local frequencies.

**Output**
- gsm.csv: contains all selected GSMs and their frequencies in all reference genomes.
- genomeGSM: a directory contains all selected GSMs and frequencies, splited into each genome. 

### Counting GSM frequencies in reads

```
go run kmers_read.go readfile GSMfile GSM_dir K
```
**Input**
- readfile: file contains reads.
- GSMfile: gsm.csv is an output file from **kmers_selection.go**.
- GSM_dir: genomeGSM is an output directory from **kmers_selection.go**.
- K: length of marker.

**Output**
- b: a directory contains frequencies of selected GSM in reads, splited into each genome where they occur.
