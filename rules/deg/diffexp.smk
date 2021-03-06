#note conda and script paths are relative paths to the rule directory
# e.g viper/rules/deg/diffexp.smk to viper/scripts/FeatCountSummary.R  is referred by ../../scripts/FeatCountSummary.R


############## deseq2

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


# rule deseq2_init:
#     input:
#         counts="results/quantification/counts/{ref}/allFC.rds"#"results/quantification/counts/all.tsv"
#     output:
#         "results/deg/deseq2/{ref}/all.rds"
#     params:
#         samples=config["samples"],
#         units=config["units"]
#     conda:
#         "../../envs/r35.yaml"
#     log:
#         "logs/deseq2/{ref}/init.log"
#     threads: get_deseq2_threads()
#     script:
#         "../../scripts/deseq2-init.R"

# rule pca:
#     input:
#         "results/deg/deseq2/{ref}/all.rds"
#     output:
#         report("results/plot/deseq2/{ref}/pca.svg", "../../report/pca.rst")
#     params:
#         pca_labels=config["pca"]["labels"]
#     conda:
#         "../../envs/r35.yaml"
#     log:
#         "logs/deseq2/{ref}/pca.log"
#     script:
#         "../../scripts/plot-pca.R"


# rule deseq2:
#     input:
#         "results/deg/deseq2/{ref}/all.rds"
#     output:
#         table=report("results/plot/deseq2/{ref}/diffexp/{contrast}.diffexp.tsv", "../../report/diffexp.rst"),
#         ma_plot=report("results/plot/deseq2/{ref}/diffexp/{contrast}.ma-plot.svg", "../../report/ma.rst"),
#     params:
#         contrast=get_contrast
#     conda:
#         "../../envs/r35.yaml"
#     log:
#         "logs/deseq2/{ref}/{contrast}.diffexp.log"
#     threads: get_deseq2_threads
#     script:
#         "../../scripts/deseq2.R"

############## edegR

rule edegR_deg:
    input:
        # counts="results/quantification/counts/{ref}/allFC.rds"
        counts="results/quantification/counts/{ref}/{readtype}/{ctype}/{RDStype}.rds",
        viper=rules.install_R_package_viper.output
    output:
        "results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_ResData.rds",
        "results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_Res.csv",
        "results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_Res.rds",
        "results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_Res_sig_"+config["diffexp"]["sig"]+"_MYlog2FC_"+config["diffexp"]["log2FC"]+".csv",
        "results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_Res_sig_"+config["diffexp"]["sig"]+".csv"
    params:
        samples=config["samples"],
        units=config["units"],
        contrast=get_contrast,
        sig=config["diffexp"]["sig"],
        log2FC=config["diffexp"]["log2FC"],
        MinGeneReads=config["diffexp"]["MinGeneReads"] #if 0 then the filter is disabled  
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_init.log"
    threads: get_deseq2_threads()
    script:
        "../../scripts/edegR-deg.R"

rule edegR_plot_global:
    input:
        counts="results/quantification/counts/{ref}/{readtype}/{ctype}/{RDStype}.rds",
        viper=rules.install_R_package_viper.output
    output:
        report("results/plot/edegR/{ref}_{readtype}_{ctype}_{RDStype}_edegR_plot_MDS.png", caption="../../report/edegR_plot_MDS.rst", category="MDS"),
        report("results/plot/edegR/{ref}_{readtype}_{ctype}_{RDStype}_edegR_plot_PCA.png", caption="../../report/pca.rst", category="PCA"),
        report("results/plot/edegR/{ref}_{readtype}_{ctype}_{RDStype}_edegR_plot_biotype.png", caption="../../report/edegR_plot_biotype.rst", category="Pie biotype")
    params:
        samples=config["samples"],
        units=config["units"]
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/edegR/{ref}/{readtype}/{ctype}/{RDStype}_edegR_plot_global.log"
    threads: get_deseq2_threads()
    script:
        "../../scripts/edegR-plot.R"

rule edegR_ResSigFC_pair_R:
    input:
        e1="results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast1}_edegR_ResData.rds",
        e2="results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast2}_edegR_ResData.rds",
        viper=rules.install_R_package_viper.output
    output:    
        o1="results/plot/edegR/{ref}_{readtype}/{ctype}_{RDStype}_ResSiglog2FC/{contrast1}_{contrast2}/VennSet.pdf",
        # o2="results/plot/edegR/{ref}_{readtype}/{ctype}_{RDStype}_ResSiglog2FC/{contrast1}_{contrast2}/Genes_Filter__{contrast1}__{contrast2}.rds",
        o3="results/plot/edegR/{ref}_{readtype}/{ctype}_{RDStype}_ResSiglog2FC/{contrast1}_{contrast2}/Genes__{contrast1}__{contrast2}.rds"
    params:
        MinGeneReads=config["diffexp"]["MinGeneReads"],
        MeanReads=1, ### Todo move param to config
        doCorrelationAnalysis="FALSE" ### Todo move param to config
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/plot/edegR/{ref}_{readtype}_{ctype}_{RDStype}_ResSiglog2FC_{contrast1}_{contrast2}.log"
    script:
        "../../scripts/edegR-deg_res_pairwise.R"       
