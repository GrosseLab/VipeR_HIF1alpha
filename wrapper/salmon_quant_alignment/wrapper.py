__author__ = "Weinholdt Claus"
__copyright__ = "Copyright 2018, Weinholdt Claus"
__email__ = "claus.weinholdt@informatik.uni-halle.de."
__license__ = "MIT"

import os
from os import path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
zip_extension = snakemake.params.get("zip_extension", "")
libtype = snakemake.params.get("libtype", "A")

bam = snakemake.input.get("bam")
transcriptome =  snakemake.input.get("transcriptome")
gtf = snakemake.input.get("gtf")

outdir = path.dirname(snakemake.output.get('quant'))

shell("salmon quant "
	" --targets {transcriptome} "
	"--alignments {bam} " 
	"--geneMap {gtf} "
	"-l {libtype} -o {outdir} "
	"-p {snakemake.threads} {extra} {log} ")

shell("cp {outdir}/libParams/flenDist.txt {outdir}/aux_info/ ")
shell("cp {outdir}/libParams/salmon_quant.log {outdir}/aux_info/ ")
