
import subprocess

### from martin 
def getParams(program, parameterSet):
    parameterSet = parameterSet.lstrip(' -')
    return(" ".join(config["parameters"][parameterSet][program]) if (program and program in config.get("parameters", {}).get(parameterSet, {}))  else "")

wildcard_constraints:
  parameterSet="(-{1,2}.+)?",
  parameterSet2="(-{1,2}.+)?",
  mapper="[^- \/]+",
  ref="[^\/]+",
  kind="[^\/]+",
  sample="[^\/]+"

rule gzip:
    """gzip any file for compression"""
    input: "{file}"
    output: "{file}.gz"
    threads: 20
    shell: "pigz -f -c -p {threads} {input} > {output}"

rule build_R_package_viper:
    input: "viper" #"{pack}"
    output: "viper.build.txt" # "{pack}_{vers}.tar.gz"
    conda: "../envs/r35.yaml"
    log: "logs/build_R_package_viper.log"#"logs/build_R_package_{pack}_{vers}.log"    
    priority: 100
    shell: 
        "R CMD build {input} && "
        " touch {output} "
        # "R CMD check {output} --no-manual --no-build-vignettes "

rule install_R_package_viper:
    input: rules.build_R_package_viper.output #"{pack}_{vers}.tar.gz"
    output: "viper.txt" #"{pack}_{vers}.txt"
    conda: "../envs/r35.yaml"
    log: "logs/build_R_install_viper.log" #"logs/build_R_install_{pack}_{vers}.log"    
    priority: 100
    script:
        "../scripts/installPack.R"


##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


# def get_trimmed_reads(wildcards):
#     """Get trimmed reads of given sample-unit."""
#     if not is_single_end(**wildcards):
#         # paired-end sample
#         return expand("results/trimmed/cutadapt/{sample}-{unit}.{group}.fastq.gz",
#                       group=[1, 2], **wildcards)
#     # single end sample
#     return "results/trimmed/cutadapt/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_trimmed_cutadapt(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        # paired-end sample
        return expand("results/trimmed/cutadapt/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "results/trimmed/cutadapt/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_trimmed_sickle(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        # paired-end sample
        return expand("results/trimmed/sickle/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "results/trimmed/sickle/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)

def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]    