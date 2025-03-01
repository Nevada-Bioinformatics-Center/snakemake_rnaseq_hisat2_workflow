import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
adapter=config["ref"]["adapter"]
#print(units)


##### target rules #####

rule all:
    input:
        #expand("logs/trimmomatic/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        directory("index_genome"),
        expand("qc/fastqc_pretrim/{unit.sample}-{unit.unit}_r1.html", unit=units.itertuples()),
        expand("qc/fastqc_pretrim/{unit.sample}-{unit.unit}_r2.html", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/{unit.sample}-{unit.unit}_r1.html", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/{unit.sample}-{unit.unit}_r2.html", unit=units.itertuples()),
        "qc/multiqc_report_pretrim.html",
        "qc/multiqc_report_posttrim.html",
        "qc/multiqc_report_all.html",
        #expand("mapped/{unit.sample}-{unit.unit}.sorted.bam", unit=units.itertuples()),
        "results/featureCounts/all.featureCounts",
        #"results/featureCounts/all_multi.featureCounts",
        #"results/featureCounts/all_multi_frac.featureCounts"

##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/fpkm.smk"
