rule david_ResSig:
    input:
        i="results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_Res_sig_"+config["diffexp"]["sig"]+".csv"
    output:
        # o=directory("results/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast}_{davidVersion}"),
        o="results/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast}_{davidVersion}_ResSig" #/{davidVersion}_chartReport.txt"
    params:
        "-t ENSEMBL_GENE_ID"
        # lambda wildcards: "-t {}".format(config["ref"][wildcards.ref]["annotation"])
    conda:
        "../../envs/py2.yaml"
    log:
        "logs/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast}_{davidVersion}_ResSig.log"
    shell:
        "python2 ./viper/scripts/DAVIDnewPythonClient/DAVID.py -i {input.i} -o {output.o} -d {wildcards.davidVersion} {params} "  

rule david_ResSigFC:
    input:
        i="results/deg/edegR/{ref}/{readtype}/{ctype}/{RDStype}_{contrast}_edegR_Res_sig_"+config["diffexp"]["sig"]+"_MYlog2FC_"+config["diffexp"]["log2FC"]+".csv"
    output:
        # o=directory("results/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast}_{davidVersion}_ResSiglog2FC")#/{davidVersion}_chartReport.txt"
        o="results/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast}_{davidVersion}_ResSiglog2FC/{davidVersion}_chartReport.txt",
        o2="results/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast}_{davidVersion}_ResSiglog2FC/{davidVersion}_chartReport_T1.txt"
    params:
        "-t ENSEMBL_GENE_ID"
        # lambda wildcards: "-t {}".format(config["ref"][wildcards.ref]["annotation"])
    conda:
        "../../envs/py2.yaml"
    log:
        "logs/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast}_{davidVersion}_ResSiglog2FC.log"
    shell:
        "python2 ./viper/scripts/DAVIDnewPythonClient/DAVID.py -i {input.i} -o {output.o} -d {wildcards.davidVersion} {params} "  


rule david_ResSigFC_pair_R:
    input:
        TrGe="results/quantification/counts/{ref}/TrGe.rds",
        D1="results/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast1}_{davidVersion}_ResSiglog2FC/{davidVersion}_chartReport_T1.txt",
        D2="results/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{contrast2}_{davidVersion}_ResSiglog2FC/{davidVersion}_chartReport_T1.txt",
        viper=rules.install_R_package_viper.output
    output:    
        o1="results/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{davidVersion}_ResSiglog2FC/{contrast1}_{contrast2}/ChartReportSigList.RDS"
    params:
        ### not needed. can get by R directly from config 
        # samples=config["samples"],
        # units=config["units"],
        # contrastNames= lambda wildcards: "{wildcards.contrast1},{wildcards.contrast2}",
        # sig=config["diffexp"]["sig"],
        # log2FC=config["diffexp"]["log2FC"]
        ""
    conda:
        "../../envs/r35.yaml"
    log:
        "logs/annotation/DAVID/{ref}_{readtype}_{ctype}_{RDStype}_{davidVersion}_ResSiglog2FC_{contrast1}_{contrast2}.log"
    script:
        "../../scripts/DAVID_output_analysis.R"       