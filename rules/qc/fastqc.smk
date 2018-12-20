def get_fastq_fq1(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()

def get_fastq_fq2(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()

rule fastqc:
    input:
        get_fastq
    output:
        html="results/qc/fastqc/{sample}-{unit}.html",
        zip="results/qc/fastqc/{sample}-{unit}_fastqc.zip"
    wrapper:
        "file:viper/wrapper/fastqc_v0.30.0"
        # "0.30.0/bio/fastqc"

rule fastqc_q1:
    input:
        get_fastq_fq1
    output:
        html="results/qc/fastqc/{sample}-{unit}_R1.html",
        zip="results/qc/fastqc/{sample}-{unit}_R1.zip"
    params: "-k 10"
    threads: 10
    log:
        "logs/fastqc/{sample}-{unit}_R1.log"
    wrapper:
        "file:viper/wrapper/fastqc_v0.17.4" ###add threads to fastqc
        # "file:viper/wrapper/fastqc_v0.30.0"
        # "0.30.0/bio/fastqc"

rule fastqc_q2:
    input:
        get_fastq_fq2
    output:
        html="results/qc/fastqc/{sample}-{unit}_R2.html",
        zip="results/qc/fastqc/{sample}-{unit}_R2.zip"
    params: "-k 10"
    threads: 10
    log:
        "logs/fastqc/{sample}-{unit}_R2.log"
    wrapper:
        "file:viper/wrapper/fastqc_v0.17.4" ###add threads to fastqc
        # "0.30.0/bio/fastqc"

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

