rule star_index:
	input:
		fasta = lambda wildcards: config["ref"][wildcards.ref]["genome"],
		gtf = lambda wildcards: config["ref"][wildcards.ref]["annotation"]
	output:
		dir = directory("references/{ref}/STAR_INDEX")
	log:
		"logs/STAR/{ref}/star_index.log"
	threads: 20
	params:
		# optional parameters
		sjdbGTFfeatureExon = lambda wildcards: config["ref"][wildcards.ref]["sjdbGTFfeatureExon"],
		sjdbGTFtagExonParentTranscript = lambda wildcards: config["ref"][wildcards.ref]["sjdbGTFtagExonParentTranscript"],
		sjdbGTFtagExonParentGene = lambda wildcards: config["ref"][wildcards.ref]["sjdbGTFtagExonParentGene"],
		sjdbOverhang = lambda wildcards: config["ref"][wildcards.ref]["sjdbOverhang"],
		limitGenomeGenerateRAM = '80000000000',
		extra=""
	wrapper:
		"file:viper/wrapper/star_index"

#/home/adsvy/anaconda3/bin/snakemake --use-conda --core 15 -kp
#/home/adsvy/anaconda3/bin/snakemake --use-conda -n
