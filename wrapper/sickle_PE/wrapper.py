__author__ = "Weinholdt Claus"
__copyright__ = "Copyright 2018, Weinholdt Claus"
__email__ = "claus.weinholdt@informatik.uni-halle.de."
__license__ = "MIT"

import os
from snakemake.shell import shell
from os import path

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

r1 = snakemake.input.get("r1")
r2 = snakemake.input.get("r2")

outdir = path.dirname(snakemake.output.get('r1'))

shell("pigz -dc -p {snakemake.threads} {r1} > {snakemake.output.r1tmp}")
shell("pigz -dc -p {snakemake.threads} {r2} > {snakemake.output.r2tmp}")

shell(
	"(sickle pe -f {snakemake.output.r1tmp} -r {snakemake.output.r2tmp} "
	"-o {snakemake.output.r1} -p {snakemake.output.r2} "
	"-s {snakemake.output.rs} -t {snakemake.params.qual_type} "
	"{extra}) {log}"
)

# shell("pigz -f -p {snakemake.threads} {snakemake.output.r1}")
# shell("pigz -f -p {snakemake.threads} {snakemake.output.r2}")
# shell("pigz -f -p {snakemake.threads} {snakemake.output.rs}")



################################################################
# def manual_decompression (reads, zip_ext):
# 	""" Allow *.bz2 input into salmon. Also provide same
# 	decompression for *gz files, as salmon devs mention
# 	it may be faster in some cases."""
# 	if zip_ext and reads:
# 		if zip_ext == 'bz2':
# 			reads = ' < (bunzip2 -c ' + reads + ')'
# 		elif zip_ext == 'gz':
# 			reads = ' < (gunzip -c ' + reads + ')'
# 	return reads

# assert (r1 is not None and r2 is not None) or r is not None, "either r1 and r2 (paired), or r (unpaired) are required as input"
# if r1:
# 	assert len(r1) == len(r2), "input-> equal number of files required for r1 and r2"
# 	if r1[0].endswith(".gz"):
# 		zip_extension = "gz"
# 		r1_cmd = ' -f ' + manual_decompression(" ".join(r1), zip_extension)
# 		r2_cmd = ' -r ' + manual_decompression(" ".join(r2), zip_extension)
# 	else:
# 		zip_extension = ""
# 		r1_cmd = ' -f ' + r1
# 		r2_cmd = ' -r ' + r2
	
# 	read_cmd = " ".join([r1_cmd,r2_cmd])
# 	shell(
# 		"(sickle pe {read_cmd} "
# 		"-o {snakemake.output.r1} -p {snakemake.output.r2} "
# 		"-s {snakemake.output.rs} -t {snakemake.params.qual_type} "
# 		"{extra}) {log}"
# 	)

# 	if r1[0].endswith(".gz")
# 	zip_extension = "gz"
# 		shell("pigz -c -p {snakemake.threads} {snakemake.output.r1} > {snakemake.output.r1}.gz")
# 		shell("pigz -c -p {snakemake.threads} {snakemake.output.r2} > {snakemake.output.r2}.gz")
# 		shell("pigz -c -p {snakemake.threads} {snakemake.output.rs} > {snakemake.output.rs}.gz")

# # if r:
# #     assert r1 is None and r2 is None, "Salmon cannot quantify mixed paired/unpaired input files. Please input either r1,r2 (paired) or r (unpaired)"
# #     r = [snakemake.input.r] if isinstance(snakemake.input.r, str) else snakemake.input.r
# #     read_cmd = ' -r ' + manual_decompression(" ".join(r), zip_extension)






