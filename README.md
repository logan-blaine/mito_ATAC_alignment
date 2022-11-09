# Alignment + preprocessing workflow for SHARE-seq data

This project uses Snakemake to align and extract mitochondrial reads + cell barcodes from raw scATAC/SHARE-seq read data whose naming conventions match the ones from https://github.com/masai1116/SHARE-seq-alignmentV2. Compatible with downstream use of [mgatk](https://github.com/caleblareau/mgatk).

## Usage examples

Run pipeline locally with 16 cores
```
snakemake --cores 16
```

## Inputs
- A folder of `.fastq` files containing paired-end reads from each sequencing run (folder name = `raw`) 
    - For each sample, files must be named "raw/{sample}.R1.fastq.gz", "raw/{sample}.R2.fastq.gz"
    - Cell barcodes extracted from each read name, currently everything after the '_'
- A file containing a list of sample names, one sample per line (e.g. `samples.txt`)
- A reference genome indexed with `bwa index`
    - Preferably with nuclear regions w/ mito homology [blacklisted](https://github.com/caleblareau/mitoblacklist)
- A config file `config.yaml` pointing to the sample file and reference genome (see example in repo)

This project uses Snakemake to perform alignment + preprocessing of scATAC (or [SHARE-seq](https://github.com/masai1116/SHARE-seq-alignmentV2/)) data, ideally for use with [my off-label mgatk fork](https://github.com/logan-blaine/mgatk/tree/barcode_fix)

## Outputs
- Aligned `.bam` files processed as follows:
    - Reads aligning to chrM only
    - Barcodes extracted to the 'CB' tag
    - Reads by read position and indexed
    - Duplicate reads removed (same barcode and 5' start positions)
    
## Requirements
Version numbers below are recommended.
```
samtools >= 1.16.1
snakemake >= 7.17.1
bwa-mem >= 0.7.17
awk (gawk) >= 5.0.1
```
