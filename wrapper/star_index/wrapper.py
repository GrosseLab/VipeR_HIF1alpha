__author__ = "Weinholdt Claus"
__copyright__ = "Copyright 2018, Weinholdt Claus"
__email__ = "claus.weinholdt@informatik.uni-halle.de."
__license__ = "MIT"

import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

print(snakemake.output.dir)
print(snakemake.params.sjdbGTFfeatureExon)
print(snakemake.params.sjdbGTFtagExonParentTranscript)
print(snakemake.params.sjdbGTFtagExonParentGene)
print(snakemake.input.gtf)
print(snakemake.params.sjdbOverhang)
print(snakemake.input.fasta)


shell("rm -fr {snakemake.output.dir}")
		
shell("mkdir {snakemake.output.dir}")    

shell(
	"STAR "
	"--runMode genomeGenerate "
	"--genomeDir {snakemake.output.dir} "
	"--genomeFastaFiles {snakemake.input.fasta} "
	"--sjdbGTFfile {snakemake.input.gtf} "
	"--sjdbGTFfeatureExon {snakemake.params.sjdbGTFfeatureExon} "
	"--sjdbGTFtagExonParentTranscript {snakemake.params.sjdbGTFtagExonParentTranscript} "
	"--sjdbGTFtagExonParentGene {snakemake.params.sjdbGTFtagExonParentGene} "
	"--sjdbOverhang {snakemake.params.sjdbOverhang} "
	"--limitGenomeGenerateRAM {snakemake.params.limitGenomeGenerateRAM} "
	"--runThreadN {snakemake.threads} "
	"{extra} "
	"{log}")

# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
# --sjdbGTFfeatureExon
# default: exon
# string: feature type in GTF file to be used as exons for building transcripts
# --sjdbGTFtagExonParentTranscript
# default: transcript id
# string: tag name to be used as exons’ transcript-parents (default
# ”transcript id” works for GTF files)
# --sjdbGTFtagExonParentGene
# default: gene id
# string: tag name to be used as exons’ gene-parents (default ”gene id” works for
# GTF files)
# --sjdbOverhang
# default: 100
# int>0: length of the donor/acceptor sequence on each side of the junctions,
# ideally = (mate length - 1)

