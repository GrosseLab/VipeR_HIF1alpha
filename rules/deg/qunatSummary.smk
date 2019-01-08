##source activate r35-ENVc3
##conda env export > r35c3.yaml

############## summarize counts 

# rule count_matrix:
#     input:
#         expand("results/mapping/star/{{ref}}/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
#     output:
#         "results/quantification/counts/{ref}/all.tsv"
#     params:
#         units=units
#     conda:
#         "../../envs/pandas.yaml"
#     script:
#         "../../scripts/count-matrix.py"

rule FeatCountSummary:
    input:
        expand("results/quantification/featureCounts/{{ref}}/{unit.sample}-{unit.unit}.feature_counts.result", unit=units.itertuples())
    output:
        "results/quantification/counts/{ref}/allFC.csv",
        "results/quantification/counts/{ref}/allFC.rds"
    params:
        units=units
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/quantification/counts/{ref}/FeatCountSummary.log"
        #"logs/deseq2/FeatCountSummary.log"
    script:
        "../../scripts/FeatCountSummary.R"


rule TximporTrGe:
    input:
        gtf = lambda wildcards: config["ref"][wildcards.ref]["annotation"],
        viper=rules.install_R_package_viper.output
    output:
        "results/quantification/counts/{ref}/TrGe.csv",
        "results/quantification/counts/{ref}/TrGe.rds"
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/quantification/counts/{ref}/TximporTrGe.log"
    script:
        "../../scripts/gtf_to_TrGe.R"

rule FeatCountSummarySetType:
    input:
        files = expand("results/quantification/featureCounts/{{ref}}/{{readtype}}/{{ctype}}/{unit.sample}-{unit.unit}/feature_counts.result", unit=units.itertuples()),
        TrGe = "results/quantification/counts/{ref}/TrGe.rds",
        viper=rules.install_R_package_viper.output
    output:
        "results/quantification/counts/{ref}/{readtype}/{ctype}/count.csv",
        "results/quantification/counts/{ref}/{readtype}/{ctype}/count.rds"
    params:
        units=units
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/quantification/counts/{ref}/{readtype}/{ctype}/FeatCountSummary.log"
        #"logs/deseq2/FeatCountSummary.log"
    script:
        "../../scripts/FeatCountSummary.R"

rule TximportData:
    input:
        files = expand("results/quantification/{{SalType}}/{{ref}}/{unit.sample}-{unit.unit}/quant.sf", unit=units.itertuples()),
        TrGe = "results/quantification/counts/{ref}/TrGe.rds",
        viper=rules.install_R_package_viper.output
    output:
        "results/quantification/counts/{ref}/{readtype}/{SalType}/estcount.rds" ### readtype currently not used for salmon
    params:
        units=units
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/quantification/counts/{ref}/{readtype}/{SalType}/TximportData.log"
    script:
        "../../scripts/TxDataSummary.R"
