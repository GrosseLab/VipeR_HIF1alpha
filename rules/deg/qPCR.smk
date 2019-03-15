rule qPCR_Anova_Normoxie:
    input:
        qPCR=config["qPCR"],
        viper=rules.install_R_package_viper.output
    output:
        "results/plot/qPCR/Normoxia/Anova_Normoxie_Q-GLC.pdf",
        "results/qPCR/Normoxia/Anova_Normoxie_Q-GLC_res.csv",
        "results/qPCR/Normoxia/Anova_Normoxie_Q-GLC_res.rds",
        report("results/plot/qPCR/Normoxia/Anova_Normoxie_Q-GLC_LDHA_ANOVA.png", caption="../../report/qPCR_res_LDHA.rst", category="pPCR")
    message:
            "`qPCR_analysis.RDStype` is project specific"
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/qPCR/Anova_Normoxie_Q-GLC.log"
    script:
        "../../scripts/qPCR_analysis.R"


rule qPCR_Normoxia_Hypoxia:
    input:
        qPCR2=config["qPCR2"],
        viper=rules.install_R_package_viper.output
    output:
        png=report("results/plot/qPCR/Normoxia_Hypoxia/LDHA.png", caption="../../report/qPCR2_res_LDHA.rst", category="pPCR")
    message:
        "`qPCR_data2_analysis.R` is project specific"    
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/qPCR/Normoxia_Hypoxia.log"
    script:
        "../../scripts/qPCR_data2_analysis.R"

    

rule qPCR_Normoxia_Hypoxia_RNAseq:
    input:
        qPCR2=config["qPCR2"],
        rds1="results/plot/edegR/{ref}_{readtype}/{ctype}_{RDStype}_ResSiglog2FC/{contrast1}_{contrast2}/Genes_Filter__{contrast1}__{contrast2}.rds",
        rds2="results/plot/edegR/{ref}_{readtype}/{ctype}_{RDStype}_ResSiglog2FC/{contrast1}_{contrast2}/Genes__{contrast1}__{contrast2}.rds",
        viper=rules.install_R_package_viper.output
    output:
        png1=report("results/plot/edegR/{ref}_{readtype}/{ctype}_{RDStype}_ResSiglog2FC/{contrast1}_{contrast2}/Genes_Filter__{contrast1}__{contrast2}_qPCR_scatterplot.png", caption="../../report/qPCR_RNAseq.rst", category="pPCR"),
        png2="results/plot/edegR/{ref}_{readtype}/{ctype}_{RDStype}_ResSiglog2FC/{contrast1}_{contrast2}/Genes__{contrast1}__{contrast2}_qPCR_scatterplot.png"
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/qPCR/Normoxia_Hypoxia_{ref}_{readtype}_{ctype}_{RDStype}_ResSiglog2FC_{contrast1}_{contrast2}_scatterplot.log"
    message: 
        "`qPCR_data2_to_RNAseq_analysis.R` is project specific" 
    script:
        "../../scripts/qPCR_data2_to_RNAseq_analysis.R"