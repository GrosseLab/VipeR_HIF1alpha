#note conda and script paths are relative paths to the rule directory
# e.g viper/rules/deg/diffexp.smk to viper/scripts/FeatCountSummary.R  is referred by ../../scripts/FeatCountSummary.R


############## deseq2

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        counts="results/quantification/counts/{ref}/allFC.rds"#"results/quantification/counts/all.tsv"
    output:
        "results/deg/deseq2/{ref}/all.rds"
    params:
        samples=config["samples"],
        units=config["units"]
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/deseq2/{ref}/init.log"
    threads: get_deseq2_threads()
    script:
        "../../scripts/deseq2-init.R"

rule pca:
    input:
        "results/deg/deseq2/{ref}/all.rds"
    output:
        report("results/plot/deseq2/{ref}/pca.svg", "../../report/pca.rst")
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/deseq2/{ref}/pca.log"
    script:
        "../../scripts/plot-pca.R"


rule deseq2:
    input:
        "results/deg/deseq2/{ref}/all.rds"
    output:
        table=report("results/plot/deseq2/{ref}/diffexp/{contrast}.diffexp.tsv", "../../report/diffexp.rst"),
        ma_plot=report("results/plot/deseq2/{ref}/diffexp/{contrast}.ma-plot.svg", "../../report/ma.rst"),
    params:
        contrast=get_contrast
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/deseq2/{ref}/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../../scripts/deseq2.R"

############## edegR

rule edegR_deg:
    input:
        # counts="results/quantification/counts/{ref}/allFC.rds"
        counts="results/quantification/counts/{ref}/{readtype}/{ctype}/{RDStype}.rds"
    output:
        "results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_ResData.rds",
        "results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_Res.csv",
        "results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_Res.rds"
    params:
        samples=config["samples"],
        units=config["units"],
        contrast=get_contrast
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_init.log"
    threads: get_deseq2_threads()
    script:
        "../../scripts/edegR-deg.R"


