rule make_transciptome_exon_fasta:
	input:
		gtf = lambda wildcards: config["ref"][wildcards.ref]["annotation"],
		genome = lambda wildcards: config["ref"][wildcards.ref]["genome"] #,viper=rules.install_R_package_viper.output
	output:
		outputFile= "references/{ref}/GTF_EXON.fa" #config["ref"][wildcards.ref]["transcriptome"]
	log:
		"logs/salmon/{ref}/build_transcriptome_fasta.log"
	threads: 15
	params:
		# optional parameters
		extra=""
	conda:
		"../../envs/r35salmon.yaml"
		#"../../envs/r35.yaml"
	script:
		"../../scripts/convert_gtf_fasta.R"

rule salmon_index:
	input:
		"references/{ref}/GTF_EXON.fa" #lambda wildcards: config["ref"][wildcards.ref]["transcriptome"]
	output:
		directory("references/{ref}/salmon_transcriptome_index")
	log:
		"logs/salmon/{ref}/transcriptome_index.log"
	threads: 24
	params:
		# optional parameters
		extra=""
	wrapper:
		"file:viper/wrapper/salmon_index"
		#"0.30.0/bio/salmon/index"

rule salmon_quant_reads:
	input:
		# If you have multiple fastq files for a single sample (e.g. technical replicates)
		# use a list for r1 and r2.
		# r1 = "results/trimmed/sickle/{sample}-{unit}.1.fastq.gz",
		# r2 = "results/trimmed/sickle/{sample}-{unit}.2.fastq.gz",
		sample = get_trimmed_sickle,
		index = "references/{ref}/salmon_transcriptome_index"
	output: 
		# auxDIR = directory('results/quantification/salmonReads/{ref}/{sample}-{unit}/aux_info'),
		quant = 'results/quantification/salmonReads/{ref}/{sample}-{unit}/quant.sf',
		lib = 'results/quantification/salmonReads/{ref}/{sample}-{unit}/lib_format_counts.json'
	log:
		'logs/salmonReads/{ref}/{sample}-{unit}.log'
	params:
		# optional parameters
		libtype ="A",
		#zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
		extra="--seqBias --gcBias --dumpEq --dumpEqWeights --useVBOpt --numBootstraps 100 --writeOrphanLinks --writeUnmappedNames"
	threads: 24
	wrapper:
		"file:viper/wrapper/salmon_quant_reads"
		#"0.30.0/bio/salmon/index"


rule salmon_quant_alignment:
	input:
		bam = "results/mapping/star/{ref}/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
		transcriptome = "references/{ref}/GTF_EXON.fa", #lambda wildcards: config["ref"][wildcards.ref]["transcriptome"],
		gtf = lambda wildcards: config["ref"][wildcards.ref]["annotation"]
	output:
		# auxDIR = directory('results/quantification/salmonAlignment/{ref}/{sample}-{unit}/aux_info'),
		quant = "results/quantification/salmonAlignment/{ref}/{sample}-{unit}/quant.sf"
	log:
		'logs/salmonAlignment/{ref}/{sample}-{unit}.log'
	threads: 24
	params:
		# optional parameters
		libtype ="A",
		extra="--seqBias --gcBias --dumpEq --dumpEqWeights --useVBOpt --numBootstraps 100 --writeOrphanLinks --writeUnmappedNames"
	wrapper:
		"file:viper/wrapper/salmon_quant_alignment"
		#"0.30.0/bio/salmon/index"