def get_bams(wildcards):
    sam = expand("star/{sample}-{unit}/Aligned.sortedByCoord.out.bam", **wildcards)
    print(sam)
    return sam

rule featurecounts_onefile:
    input:
        sam=expand("mapped/{unit.sample}.{unit.unit}.bam", unit=units.itertuples()),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/featureCounts/all",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    log:
        "logs/featurecount/all.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecounts"])
    threads: 8
    wrapper:
        "0.69.0/bio/subread/featurecounts"

rule featurecounts_onefile_multimap:
    input:
        sam=expand("mapped/{unit.sample}.{unit.unit}.bam", unit=units.itertuples()),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/featureCounts/all_multi",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    log:
        "logs/featurecount/all_multi.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecountsmulti"])
    threads: 8
    wrapper:
        "0.69.0/bio/subread/featurecounts"

rule featurecounts_onefile_multimap_frac:
    input:
        sam=expand("mapped/{unit.sample}.{unit.unit}.bam", unit=units.itertuples()),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/featureCounts/all_multi_frac",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    log:
        "logs/featurecount/all_multi.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecountsmulti"])
    threads: 8
    wrapper:
        "0.69.0/bio/subread/featurecounts"

#rule insert_size:
#    input:
#        #"star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
#        "star/{sample}-{unit}/{sample}-{unit}_Aligned.sortedByCoord.out.bam",
#    output:
#        txt="stats/{sample}-{unit}.isize.txt",
#        pdf="stats/{sample}-{unit}.isize.pdf"
#    log:
#        "logs/picard/insert_size/{sample}-{unit}.log"
#    params:
#        # optional parameters (e.g. relax checks as below)
#        "VALIDATION_STRINGENCY=LENIENT "
#        #"METRIC_ACCUMULATION_LEVEL=null "
#        #"METRIC_ACCUMULATION_LEVEL=SAMPLE"
#    # optional specification of memory usage of the JVM that snakemake will respect with global
#    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
#    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
#    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
#    resources:
#        mem_mb=1024
#    wrapper:
#        "0.68.0/bio/picard/collectinsertsizemetrics"

