##Usage

### Merge GSM frequencies in genomes and reads

```
python merge_Fc.py genomeGSM b output_dir
```
**Input**
- genomeGSM: directory stores csv files which contains GSM frequencies in reference genomes. Each csv file has two columns: GSM_id and frequency.
- b: directory stores csv files which containes GSM frequencies in reads. Each csv file has two columns: GSM_id and frequency.

**Output**
- output_dir: directory stores csv files which contains GSM frequencies in reference genomes and reads. Each csv file has three columns: GSM_id, frequency in reference genome, frequency in reads.

### Estimate Absolute Abundance

```
python le_abundance.py Fc_dir numberOfRefGenomes output_file
```

**Input**
- Fc_dir: directory stores csv files which contains GSM frequencies in reference genomes and reads. Output directory from **merge_Fc.py**
- numberOfRefGenomes: number of reference genomes in database.

**Output**
- output_file: contains the estimated absolute abundances, should be in csv format. This file has two columns: genome name and its estimated absolute abundance. 
