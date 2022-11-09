# Alignment + preprocessing workflow for SHARE-seq data

This project uses Snakemake to align mitochondrial reads and extract cell barcodes from scATAC/SHARE-seq read data processed using [Sai Ma's pipeline](https://github.com/masai1116/SHARE-seq-alignmentV2). Compatible with downstream use of [my mgatk fork](https://github.com/logan-blaine/mgatk).

## Usage examples

Run pipeline locally with 16 cores (alignment workflow only)
```
snakemake --cores 16
```

Run pipeline locally with 16 cores (alignment workflow + mgatk)
```
snakemake --cores 16 all
```

## Inputs
- A folder of `.fastq` files containing paired-end reads from each sequencing run (folder name = `raw`) 
    - For each sample, files must be named "raw/{sample}.R1.fastq.gz", "raw/{sample}.R2.fastq.gz"
    - Cell barcodes extracted from each read name, currently everything after the '_'
- A file containing a list of sample names, one sample per line (e.g. `samples.txt`)
- A reference genome indexed with `bwa index`
    - Preferably with nuclear regions w/ mito homology [blacklisted](https://github.com/caleblareau/mitoblacklist)
- A config file `config.yaml` pointing to the sample file and reference genome (see example in repo)

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
