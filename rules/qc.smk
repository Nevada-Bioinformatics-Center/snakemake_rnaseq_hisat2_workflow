def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_fastq1(wildcards):
    fq1 = units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    return fq1

def get_fastq2(wildcards):
    fq2 = units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()
    return fq2

## RSEQC


rule rseqc_gtf2bed:
    input:
        config["ref"]["annotationrseqc"],
    output:
        bed="qc/rseqc/annotation.bed",
        db=temp("qc/rseqc/annotation.db"),
    log:
        "logs/rseqc_gtf2bed.log",
    resources: time_min=120, mem_mb=8000, cpus=1
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


#rule rseqc_junction_annotation:
#    input:
#        bam="mapped/{sample}.{unit}.bam",
#        bed="qc/rseqc/annotation.bed",
#    output:
#        "qc/rseqc/{sample}-{unit}.junctionanno.junction.bed",
#    priority: 1
#    log:
#        "logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.log",
#    params:
#        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
#        prefix=lambda w, output: strip_suffix(output[0], ".junction.bed"),
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
#        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="mapped/{sample}.{unit}.sorted.bam",
        bed="qc/rseqc/annotation.bed",
    output:
        "qc/rseqc/{sample}-{unit}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.log",
    resources: time_min=120, mem_mb=8000, cpus=1
    params:
        extra=r"-q 255",
        prefix=lambda w, output: strip_suffix(
            output[0], ".junctionSaturation_plot.pdf"
        ),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "mapped/{sample}.{unit}.sorted.bam",
    output:
        "qc/rseqc/{sample}-{unit}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}-{unit}.log",
    resources: time_min=120, mem_mb=8000, cpus=1
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="mapped/{sample}.{unit}.sorted.bam",
        bed="qc/rseqc/annotation.bed",
    output:
        "qc/rseqc/{sample}-{unit}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}-{unit}.log",
    resources: time_min=120, mem_mb=8000, cpus=1
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="mapped/{sample}.{unit}.sorted.bam",
        bed="qc/rseqc/annotation.bed",
    output:
        "qc/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".inner_distance.txt"),
    resources: time_min=120, mem_mb=8000, cpus=1
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="mapped/{sample}.{unit}.sorted.bam",
        bed="qc/rseqc/annotation.bed",
        bai="mapped/{sample}.{unit}.sorted.bam.bai",
    output:
        "qc/rseqc/{sample}-{unit}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}-{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    resources: time_min=120, mem_mb=8000, cpus=1
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        bam="mapped/{sample}.{unit}.sorted.bam",
        bai="mapped/{sample}.{unit}.sorted.bam.bai"
    output:
        "qc/rseqc/{sample}-{unit}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".DupRate_plot.pdf"),
    resources: time_min=120, mem_mb=16000, cpus=1
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        bam="mapped/{sample}.{unit}.sorted.bam",
        bai="mapped/{sample}.{unit}.sorted.bam.bai"
    output:
        "qc/rseqc/{sample}-{unit}.readgc.GC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".GC_plot.pdf"),
    resources: time_min=120, mem_mb=8000, cpus=1
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input.bam} -o {params.prefix} > {log} 2>&1"

rule fastqc_pretrim_r1:
    input:
       get_fastq1
    output:
        html="qc/fastqc_pretrim/{sample}-{unit}_r1.html",
        zip="qc/fastqc_pretrim/{sample}-{unit}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}-{unit}_r1.log"
    threads: 1
    wrapper:
        "v0.69.0/bio/fastqc"

rule fastqc_pretrim_r2:
    input:
       get_fastq2
    output:
        html="qc/fastqc_pretrim/{sample}-{unit}_r2.html",
        zip="qc/fastqc_pretrim/{sample}-{unit}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}-{unit}_r2.log"
    threads: 1
    wrapper:
        "v0.69.0/bio/fastqc"

rule fastqc_posttrim_r1:
    input:
        "trimmed/{sample}-{unit}.1.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{sample}-{unit}_r1.html",
        zip="qc/fastqc_posttrim/{sample}-{unit}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{sample}-{unit}_r1.log"
    threads: 1
    wrapper:
        "v0.69.0/bio/fastqc"

rule fastqc_posttrim_r2:
    input:
        "trimmed/{sample}-{unit}.2.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{sample}-{unit}_r2.html",
        zip="qc/fastqc_posttrim/{sample}-{unit}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{sample}-{unit}_r2.log"
    threads: 1
    wrapper:
        "v0.69.0/bio/fastqc"

rule multiqc_pre:
    input:
        expand("qc/fastqc_pretrim/{unit.sample}-{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_pretrim/{unit.sample}-{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_pretrim.html"
    log:
        "logs/multiqc_pre.log"
    wrapper:
        "0.62.0/bio/multiqc"

rule multiqc_post:
    input:
        expand("logs/trimmomatic/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        #expand("report/pe/{unit.sample}-{unit.unit}.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/{unit.sample}-{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/{unit.sample}-{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_posttrim.html"
    log:
        "logs/multiqc_post.log"
    wrapper:
        "0.62.0/bio/multiqc"

rule multiqc:
    input:
        expand("mapped/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples()),
        expand("logs/hisat2/{unit.sample}_{unit.unit}.log", unit=units.itertuples()),
        expand("results/featureCounts/all.featureCounts.summary", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        expand("qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("qc/rseqc/{unit.sample}-{unit.unit}.stats.txt", unit=units.itertuples()),
        expand("qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        expand("qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("qc/rseqc/{unit.sample}-{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("qc/rseqc/{unit.sample}-{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        #expand("logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        expand("logs/trimmomatic/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        #expand("report/pe/{unit.sample}-{unit.unit}.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/{unit.sample}-{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/{unit.sample}-{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_all.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.62.0/bio/multiqc"
