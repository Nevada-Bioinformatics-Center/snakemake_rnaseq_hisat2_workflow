def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_fastq1(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna().item()

def get_fastq2(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna().item()

rule trimmomatic_pe:
    input:
        r1=get_fastq1,
        r2=get_fastq2
    output:
        r1="trimmed/{sample}-{unit}.1.fastq.gz",
        r2="trimmed/{sample}-{unit}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed/{sample}-{unit}.1.unpaired.fastq.gz",
        r2_unpaired="trimmed/{sample}-{unit}.2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    params:
        # list of trimmers (see manual)
        #trimmer = ["ILLUMINACLIP:/data/projects/common_references/raw/CustomBlacklist.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36"],
        trimmer = [f"ILLUMINACLIP:{config['ref']['adapter']}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads: 4
    resources: time_min=820, mem_mb=20000, cpus=4
    wrapper:
        "v0.69.0/bio/trimmomatic/pe"

#rule fastp_pe:
#    input:
#        sample=get_fastq
#    output:
#        trimmed=["trimmed/{sample}-{unit}.1.fastq.gz", "trimmed/{sample}-{unit}.2.fastq.gz"],
#        html="report/pe/{sample}-{unit}.html",
#        json="report/pe/{sample}-{unit}.json"
#    log:
#        "logs/fastp/pe/{sample}-{unit}.log"
#    params:
#        adapters="--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
#        extra=""
#    threads: 4
#    wrapper:
#        "0.72.0/bio/fastp"
