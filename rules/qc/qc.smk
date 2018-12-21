
rule samtools_stats:
    input:
        "results/mapping/star/{ref}/{sample}-{unit}/Aligned.out.sorted.bam"
    output:
        "results/qc/samtools-stats/{ref}/{sample}-{unit}.txt"
    log:
        "logs/samtools-stats/{ref}/{sample}-{unit}.log"
    wrapper:
        "file:viper/wrapper/samtools_v0.30.0/stats"
        #"0.30.0/bio/samtools/stats"
        
### Note multiple inputs like from featureCounts or salmon can not be distinguished -> split in independent multiqc rules
rule multiqc:
    input:
        expand(["results/trimmed/cutadapt/{u.sample}-{u.unit}.qc.txt",
                "results/qc/samtools-stats/{{ref}}/{u.sample}-{u.unit}.txt",
                # "results/qc/fastqc/{u.sample}-{u.unit}_fastqc.zip",
                "results/qc/fastqc/{u.sample}-{u.unit}_R1_fastqc.zip",
                "results/qc/fastqc/{u.sample}-{u.unit}_R2_fastqc.zip",
                "results/quantification/featureCounts/{{ref}}/PE/unique/{u.sample}-{u.unit}/feature_counts.result.summary",
                "results/mapping/star/{{ref}}/{u.sample}-{u.unit}/Log.final.out",
                "results/quantification/salmonAlignment/{{ref}}/{u.sample}-{u.unit}/aux_info"
                ],
               u=units.itertuples()),
    output:
        report("results/qc/{ref}/multiqc.html", caption="../../report/multiqc.rst", category="Quality control")
    log:
        "logs/{ref}/multiqc.log"
    wrapper:
        "file:viper/wrapper/multiqc_v0.30.0"
        # "0.30.0/bio/multiqc"

rule multiqc_salmonReads:
    input:
        expand(["results/quantification/salmonReads/{{ref}}/{u.sample}-{u.unit}/aux_info"
                ],
               u=units.itertuples()),
    output:
        report("results/qc/{ref}/multiqc_salmonReads.html", caption="../../report/multiqc.rst", category="Quality control for salmon quant Reads")
    log:
        "logs/{ref}/multiqc_salmonReads.log"
    wrapper:
        "file:viper/wrapper/multiqc_v0.30.0"
        # "0.30.0/bio/multiqc"

rule multiqc_FeatC_fraction:
    input:
        expand([
                "results/quantification/featureCounts/{{ref}}/PE/fraction/{u.sample}-{u.unit}/feature_counts.result.summary"
                ],
               u=units.itertuples()),
    output:
        report("results/qc/{ref}/multiqc_FeatC_fraction.html", caption="../../report/multiqc.rst", category="Quality control for featureCounts fraction")
    log:
        "logs/{ref}/multiqc_FeatC_fraction.log"
    wrapper:
        "file:viper/wrapper/multiqc_v0.30.0"
        # "0.30.0/bio/multiqc"

rule multiqc_FeatC_all:
    input:
        expand(["results/quantification/featureCounts/{{ref}}/PE/all/{u.sample}-{u.unit}/feature_counts.result.summary"
                ],
               u=units.itertuples()),
    output:
        report("results/qc/{ref}/multiqc_FeatC_all.html", caption="../../report/multiqc.rst", category="Quality control for featureCounts all")
    log:
        "logs/{ref}/multiqc_FeatC_all.log"
    wrapper:
        "file:viper/wrapper/multiqc_v0.30.0"
        # "0.30.0/bio/multiqc"