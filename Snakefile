configfile: "config.yaml"

sample_file = config["sample_file"]
with open(sample_file) as f:
    SAMPLES = [s.strip() for s in f.readlines()]

rule all:
    input: expand('mt-bam/{sample}.bam', sample=SAMPLES)

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
