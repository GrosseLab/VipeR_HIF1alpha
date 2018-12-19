# ruleorder: star > sambamba_sort

rule star:
	input:
		# path to STAR reference genome index
		index=lambda wildcards: directory(config["ref"][wildcards.ref]["index"]),
		sample = get_trimmed_sickle
		#fq1 = "results/trimmed/cutadapt/{sample}-{unit}.1.fastq.gz",       ### -> work with 0.30.0/bio/star/align
		#fq2 = "results/trimmed/cutadapt/{sample}-{unit}.2.fastq.gz"        ### -> work with 0.30.0/bio/star/align
	output:
		# see STAR manual for additional output files
		"results/mapping/star/{ref}/{sample}-{unit}/Aligned.out.bam",
		"results/mapping/star/{ref}/{sample}-{unit}/ReadsPerGene.out.tab",
		"results/mapping/star/{ref}/{sample}-{unit}/Aligned.toTranscriptome.out.bam"
	log:
		"logs/star/{ref}/{sample}-{unit}.log"
	message: "Executing STAR with {threads} threads on the following files {input.sample}"
	params:
		# path to STAR reference genome index
		# index=lambda wildcards: config["ref"][wildcards.ref]["index"],
		# optional parameters
		extra=lambda wildcards: "--outReadsUnmapped Fastx --outSAMattributes Standard --outSAMstrandField intronMotif --outSAMprimaryFlag AllBestScore --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile {} --sjdbGTFfeatureExon {} --sjdbGTFtagExonParentTranscript {} --sjdbOverhang {} {}".format(
			  config["ref"][wildcards.ref]["annotation"],config["ref"][wildcards.ref]["sjdbGTFfeatureExon"],config["ref"][wildcards.ref]["sjdbGTFtagExonParentTranscript"], config["ref"][wildcards.ref]["sjdbOverhang"], config["params"]["star"])
	threads: 24
	wrapper:
		"file:viper/wrapper/star_align"
		# "0.30.0/bio/star/align"
	
		### set by wrapper
		# "--outSAMtype BAM Unsorted "
		# "--outFileNamePrefix {outprefix} "
		# "--outStd Log "

# $p_ASTools/STAR-2.5.2b/bin/Linux_x86_64_static/STAR --genomeDir $GenomeDIR --readFilesIn $R1 $R2 --readFilesCommand zcat --runThreadN $Thread --outReadsUnmapped Fastx --limitGenomeGenerateRAM 800000000000  --sjdbOverhang 199 --sjdbGTFfile $gtf --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript transcript_id --alignIntronMin 20 --alignIntronMax 4000 --outFilterMismatchNmax $MaxMismatch --outFilterMultimapNmax $MaxMultimap --outSAMattributes Standard --outSAMstrandField intronMotif --outFileNamePrefix $ST_folder/ --outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode TranscriptomeSAM GeneCounts --outTmpDir /scratch/user/adsvy/_STARtmpSTAR_$genome_S.$sid > $ST_folder/accepted_hits.sort.bam
# $p_ASTools/sambamba_v0.6.3 index $ST_folder/accepted_hits.sort.bam -t 10
