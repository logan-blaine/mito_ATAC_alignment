configfile: "config.yaml"

sample_file = config["sample_file"]
with open(sample_file) as f:
    SAMPLES = [s.strip() for s in f.readlines()]

rule align:
    input: expand('mt-bam/{sample}.bam', sample=SAMPLES)
    # input: expand('mgatk-output/{sample}/final/mgatk.rds', sample=SAMPLES)
    # input: expand('mgatk-output/{sample}/final/barcodeQuants.tsv', sample=SAMPLES)

rule all: 
    input: expand('mgatk-output/{sample}/final/barcodeQuants.tsv', sample=SAMPLES)

rule bwa_mem_chrM:
    input:
        reads=["raw/{sample}.R1.fastq.gz", "raw/{sample}.R2.fastq.gz"],
        ref=config["reference"],
    output:
        "mapped/{sample}.bam",
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        awk_cmd = r'(/^@/) {print; next} ($3=="chrM") {split($1, a, "_"); printf("%s\tCB:Z:%s\n",$0,a[2])}'
    threads: 8
    shell:
        "bwa mem -t {threads} {input.ref} {input.reads} {params.extra} 2>{log}"
        " | awk \'{params.awk_cmd}\'"
        " | samtools view -b -o {output}"

rule samtools_markdup:
    input:
        "mapped/{sample}.bam",
    output:
        bam="mt-bam/{sample}.bam",
        bai="mt-bam/{sample}.bam.bai",
    params:
        tmpdir = '/tmp'
    resources:
        mem_mb = 4000
    threads: 3
    shell:
        "samtools fixmate -m {input} -"
        " | samtools sort -m {resources.mem_mb}M -T {params.tmpdir} -"
        " | samtools markdup -r --barcode-tag CB - {output.bam}"
        " && samtools index {output}"

# in progress: run mgatk on each sample 
rule mgatk_bcall:
    input:
        "mt-bam/{sample}.bam"
    output:
        # directory("mgatk-output/{sample}"), 
        "mgatk-output/{sample}/final/barcodeQuants.tsv"
    params:
        genome=config['genome'],
        outdir = lambda wildcards, output: f'mgatk-output/{wildcards["sample"]}',
        min_barcodes=config['min_barcodes'] if 'min_barcodes' in config else 100,
        barcodes=f'-b {config["barcodes"]}' if 'barcodes' in config else '',
        num_samples='-ns 800'
    threads: 8
    shell:
        "mgatk bcall -i {input} -o {params.outdir} -bt CB"
        " -mb {params.min_barcodes} {params.barcodes} {params.num_samples}"
        " -g {params.genome} -kd -c {threads}"
