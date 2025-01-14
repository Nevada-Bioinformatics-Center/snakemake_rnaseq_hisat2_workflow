def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            #print(expand("trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards))
            return expand("trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards)
        # single end sample
        #print("trimmed/{sample}-{unit}.fastq.gz".format(**wildcards))
        return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_fq(wildcards):
    if not is_single_end(**wildcards):
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards)
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_trim_fastq1(wildcards):
    fq1 = expand("trimmed/{sample}-{unit}.1.fastq.gz", **wildcards)
    #print(fq1)
    return fq1

def get_trim_fastq2(wildcards):
    fq2 = expand("trimmed/{sample}-{unit}.2.fastq.gz", **wildcards)
    #print(fq2)
    return fq2

rule hisat2_extractexons:
    input:
        gtf = config["ref"]["annotation"]
    output:
        "hisat2_prep/genome.exons"
    log:
        "logs/hisat2_extract_exons.log"
    threads: 1
    resources: time_min=120, mem_mb=2000, cpus=1
    conda:
        "../envs/hisat2.yaml"
    shell:
        "hisat2_extract_exons.py {input.gtf} > {output} 2> {log}"

rule hisat2_extractsplicesites:
    input:
        gtf = config["ref"]["annotation"]
    output:
        "hisat2_prep/genome.ss"
    log:
        "logs/hisat2_extract_ss.log"
    threads: 1
    resources: time_min=120, mem_mb=2000, cpus=1
    conda:
        "../envs/hisat2.yaml"
    shell:
        "hisat2_extract_splice_sites.py {input.gtf} > {output} 2> {log}"

rule hisat2_index:
    input:
        fasta = config["ref"]["genomefa"],
        exons = "hisat2_prep/genome.exons",
        ss = "hisat2_prep/genome.ss"
    output:
        directory("index_genome")
    params:
        prefix = "index_genome/genome"
    log:
        "logs/hisat2_index_genome.log"
    threads: 16
    resources: time_min=820, mem_mb=200000, cpus=16
    conda:
        "../envs/hisat2.yaml"
    shell:
        "mkdir {params.prefix} && hisat2-build -p {threads} --ss {input.ss} --exon {input.exons} {input.fasta} {params.prefix} 2> {log}"


rule hisat2_align:
    input:
        #reads=get_fq
        r1=get_trim_fastq1,
        r2=get_trim_fastq2,
        idx="index_genome/"
    output:
        "mapped/{sample}.{unit}.bam"
    log:
        "logs/hisat2/{sample}_{unit}.log"
    params:
      ## --new summary to allow multiqc parsing and 
      ## --dta to use XS BAM alignment information for stringtie downstream
        #extra="--new-summary --dta",
        extra="{}".format(config["params"]["hisat2"]),
        idx="index_genome/genome",
    threads: 16
    wildcard_constraints:
        unit="rep\d+"
    resources: time_min=820, mem_mb=40000, cpus=16
    conda:
        "../envs/hisat2.yaml"
    shell:
        "(hisat2 --threads {threads} -x {params.idx} {params.extra} -1 {input.r1} -2 {input.r2} | samtools view -Sbh -o {output}) 2> {log}"

rule samtools_index:
    input:
        "mapped/{sample}.{unit}.sorted.bam"
    output:
        "mapped/{sample}.{unit}.sorted.bam.bai"
    params:
        "" # optional params string
    resources: time_min=320, mem_mb=2000, cpus=1
    wrapper:
        "0.73.0/bio/samtools/index"

rule sambamba_sort:
    input:
        "mapped/{sample}.{unit}.bam"
    output:
        "mapped/{sample}.{unit}.sorted.bam"
    log:
        "logs/sambamba-sort/{sample}.{unit}.log"
    params:
        #"-m 20G" # optional additional parameters as string
        "" # optional additional parameters as string
    threads: 16 
    wildcard_constraints:
        unit="rep\d+"
    resources: time_min=1320, mem_mb=20000, cpus=16
    wrapper:
        "0.74.0/bio/sambamba/sort"

