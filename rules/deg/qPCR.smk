rule qPCR_Anova_Normoxie:
    input:
        qPCR=config["qPCR"],
        viper=rules.install_R_package_viper.output
    output:
        "results/plot/qPCR/Anova_Normoxie_Q-GLC.pdf",
        "results/qPCR/Anova_Normoxie_Q-GLC_res.csv",
        "results/qPCR/Anova_Normoxie_Q-GLC_res.rds",
        report("results/plot/qPCR/Anova_Normoxie_Q-GLC_LDHA_ANOVA.png", caption="../../report/qPCR_res_LDHA.rst", category="pPCR")
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/qPCR/Anova_Normoxie_Q-GLC.log"
    script:
        "../../scripts/qPCR_analysis.R"

