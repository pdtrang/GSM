# GSM
GSM is a tool that allows user to profile microbial genomes in the community using Genomic Specific Markers.


## Installation

1. Install Go https://golang.org/dl/
2. Install kmers package from github: http://github.com/vtphan/kmers
3. Python 2.7.x

## Reference Genomes

Reference genomes are in FASTA format. If there are more than one contigs in a genome, those contigs are concatenated and treated as a single sequence.

## Samples

User can use simulator to generate reads. Reads are in FASTA format with multiple short reads in one file. Each read is described in two lines: 
- First line: Read name/information.
- Second line: read sequence.

## Content

- main: contains Go source code for selecting GSMs and counting GSM frequencies in reads.
- solvers: contains Python source code to estimate the abundance by linear equation.


