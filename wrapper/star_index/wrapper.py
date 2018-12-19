__author__ = "Weinholdt Claus"
__copyright__ = "Copyright 2018, Weinholdt Claus"
__email__ = "claus.weinholdt@informatik.uni-halle.de."
__license__ = "MIT"

import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

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
	"--sjdbOverhang {snakemake.params.sjdbOverhang} "
	"--limitGenomeGenerateRAM {snakemake.params.limitGenomeGenerateRAM} "
	"--runThreadN {snakemake.threads} "
	"{extra} "
	"{log}")


