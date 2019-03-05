# def get_fastq(wildcards):
#     return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

rule cutadapt_pe:
	input:
		get_fastq
	output:
		fastq1="results/trimmed/cutadapt/{sample}-{unit}.1.fastq.gz",
		fastq2="results/trimmed/cutadapt/{sample}-{unit}.2.fastq.gz",
		qc="results/trimmed/cutadapt/{sample}-{unit}.qc.txt"
	params:
		"-b Illumina_Universal_Adapter=AGATCGGAAGAG "
		"-B Illumina_Universal_Adapter=AGATCGGAAGAG "
		"-b RC_Illumina_Universal_Adapter=CTCTTCCGATCT "
		"-B RC_Illumina_Universal_Adapter=CTCTTCCGATCT "
		"-b Krohn_P5=AATGATACGGCGACCACCGAGATCTACAC "
		"-b Krohn_P7=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC "
		"-B Krohn_P5=AATGATACGGCGACCACCGAGATCTACAC "
		"-B Krohn_P7=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC "
		"-a Illumina_PairedEndAdapter_1=ACACTCTTTCCCTACACGACGCTCTTCCGATCT "
		"-a Illumina_PairedEndAdapter_2=CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT "
		"-A Illumina_PairedEndAdapter_1=ACACTCTTTCCCTACACGACGCTCTTCCGATCT "
		"-A Illumina_PairedEndAdapter_2=CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT "
		"-a PrefixPE2_1=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT "
		"-A PrefixPE2_2=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT "
		"-e 0.1 -O 5 -m 40 -q 20,20 "
		"-a {} {}".format(config["adapter"] ,config["params"]["cutadapt-pe"])
	log:
		"logs/cutadapt/{sample}-{unit}.log"
	threads: 24         
	wrapper:
		"file:viper/wrapper/cutadapt_v0.30.0/pe"
		# "0.30.0/bio/cutadapt/pe"

rule cutadapt:
	input:
		get_fastq
	output:
		fastq="results/trimmed/cutadapt/{sample}-{unit}.fastq.gz",
		qc="results/trimmed/{sample}-{unit}.qc.txt"
	params:
		"-b Illumina_Universal_Adapter=AGATCGGAAGAG "
		"-b RC_Illumina_Universal_Adapter=CTCTTCCGATCT "
		"-b Krohn_P5=AATGATACGGCGACCACCGAGATCTACAC "
		"-b Krohn_P7=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC "
		"-b TruSeq_Universal_Adapter=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT "
		"-b Illumina_Multiplexing_Index_Sequencing_Primer=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
		"-b Illumina_Small_RNA_3p_Adapter_1=ATCTCGTATGCCGTCTTCTGCTTG "
		"-b Illumina_Single_End_Adapter_1=ACACTCTTTCCCTACACGACGCTGTTCCATCT "
		"-b Illumina_Single_End_Adapter_2=CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT "
		"-b Illumina_Single_End_Sequencing_Primer=ACACTCTTTCCCTACACGACGCTCTTCCGATCT "
		"-e 0.1 -O 5 -m 40 -q 20,20 "
		"-a {} {}".format(config["adapter"], config["params"]["cutadapt-se"])
	log:
		"logs/cutadapt/{sample}-{unit}.log"
	threads: 24       
	wrapper:
		"file:viper/wrapper/cutadapt_v0.30.0/se"
		# "0.30.0/bio/cutadapt/se"


ruleorder: sickle_pe > gzip
rule sickle_pe:
	input:
		r1="results/trimmed/cutadapt/{sample}-{unit}.1.fastq.gz",
		r2="results/trimmed/cutadapt/{sample}-{unit}.2.fastq.gz"
	output:
		r1="results/trimmed/sickle/{sample}-{unit}.1.fastq",
		r2="results/trimmed/sickle/{sample}-{unit}.2.fastq",
		rs="results/trimmed/sickle/{sample}-{unit}.output_single.fastq",
		# r1="results/trimmed/sickle/{sample}-{unit}.1.fastq.gz",
		# r2="results/trimmed/sickle/{sample}-{unit}.2.fastq.gz",
		# rs="results/trimmed/sickle/{sample}-{unit}.output_single.fastq.gz",
		r1tmp=temp("results/trimmed/sickle/TMP---{sample}-{unit}.1.fastq"),
		r2tmp=temp("results/trimmed/sickle/TMP---{sample}-{unit}.2.fastq")
	params:
		qual_type="sanger",
		# optional extra parameters
		extra="-q 20 -l 40 " #--gzip-output"
	log:
		# optional log file
		"results/trimmed/sickle/{sample}-{unit}.log"
	threads: 24  
	wrapper:
		"file:viper/wrapper/sickle_PE"

rule sickle_se:
	input:
		r1="results/trimmed/cutadapt/{sample}-{unit}.fastq.gz"
	output:
		r1="results/trimmed/sickle/{sample}-{unit}.fastq",
		r1tmp=temp("results/trimmed/sickle/TMP---{sample}-{unit}.fastq")
	params:
		qual_type="sanger",
		# optional extra parameters
		extra="-q 20 -l 40 "
	log:
		# optional log file
		"results/trimmed/sickle/{sample}-{unit}.log"
	wrapper:
		"file:viper/wrapper/sickle_SE"
