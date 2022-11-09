# Single-cell ATAC alignment + preprocessing for mgatk

This project uses Snakemake to perform alignment + preprocessing of scATAC (or [SHARE-seq](https://github.com/masai1116/SHARE-seq-alignmentV2/)) data, ideally for use with [my off-label mgatk fork](https://github.com/logan-blaine/mgatk/tree/barcode_fix)

## Requirements

Version numbers below are recommended.
```
samtools >= 1.16
snakemake >= 7.17.1
bwa-mem >= 0.7.17
GNU awk >= 5.0.1
```


## Inputs
- `.fastq` files containing paired-end reads from each sequencing run
    - Cell barcodes are taken from the each read name as everything after the '_'
- a reference genome
    - preferably with nuclear regions w/ mito homology [blacklisted](https://github.com/caleblareau/mitoblacklist)

## Outputs
- Aligned `.bam` files
    - Duplicates removed
    - Reads aligning to chrM only
    - Barcodes present in the 'CB' tag
