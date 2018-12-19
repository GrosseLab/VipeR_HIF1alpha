def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()

rule fastqc:
    input:
        get_fastq
    output:
        html="results/qc/fastqc/{sample}-{unit}.1.html",
        zip="results/qc/fastqc/{sample}-{unit}.1.zip"
    params: "-k 10"
    log:
        "logs/fastqc/{sample}-{unit}.1.log"
    wrapper:
        "0.30.0/bio/fastqc"

# rule cutadapt_pe:
#     input:
#         get_fastq
#     output:
#         fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
#         fastq2="trimmed/{sample}-{unit}.2.fastq.gz",
#         qc="trimmed/{sample}-{unit}.qc.txt"
#     params:
#         "-a {} {}".format(config["adapter"], config["params"]["cutadapt-pe"])
#     log:
#         "logs/cutadapt/{sample}-{unit}.log"
#     wrapper:
#         "0.17.4/bio/cutadapt/pe"


# rule cutadapt:
#     input:
#         get_fastq
#     output:
#         fastq="trimmed/{sample}-{unit}.fastq.gz",
#         qc="trimmed/{sample}-{unit}.qc.txt"
#     params:
#         "-a {} {}".format(config["adapter"], config["params"]["cutadapt-se"])
#     log:
#         "logs/cutadapt/{sample}-{unit}.log"
#     wrapper:
#         "0.17.4/bio/cutadapt/se"

