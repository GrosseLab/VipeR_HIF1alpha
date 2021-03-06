__author__ = "Weinholdt Claus"
__copyright__ = "Copyright 2018, Weinholdt Claus"
__email__ = "claus.weinholdt@informatik.uni-halle.de."
__license__ = "MIT"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

n = len(snakemake.input.sample) # just work if list is input -> for singleend we have to check if it is a string by isinstance(snakemake.input.sample, str)
print(snakemake.input.sample)
print(n)
assert isinstance(snakemake.input.sample, str) or n == 2, "input->sample must have 1 (single-end) or 2 (paired-end) elements."

if snakemake.input.sample[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
elif isinstance(snakemake.input.sample, str) and snakemake.input.sample.endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

outprefix = os.path.dirname(snakemake.output[0]) + "/"

print(readcmd)
print(outprefix)


shell(
    "STAR "
    "{extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.input.index} "      #have to set to input! so snakemake can build it.
    "--readFilesIn {snakemake.input.sample} "
    "{readcmd} "
    "--outSAMtype BAM Unsorted "
    "--outFileNamePrefix {outprefix} "
    "--outStd Log "
    "{log}")


# #https://snakemake-wrappers.readthedocs.io/en/"0.30.0/wrappers/star/align.html
#     import os
#     from snakemake.shell import shell

#     extra = snakemake.params.get("extra", "")
#     log = snakemake.log_fmt_shell(stdout=True, stderr=True)

#     fq1 = snakemake.input.get("fq1")
#     assert fq1 is not None, "input-> fq1 is a required input parameter"
#     fq1 = [snakemake.input.fq1] if isinstance(snakemake.input.fq1, str) else snakemake.input.fq1
#     fq2 =  snakemake.input.get("fq2")
#     if fq2:
#         fq2 = [snakemake.input.fq2] if isinstance(snakemake.input.fq2, str) else snakemake.input.fq2
#         assert len(fq1) == len(fq2), "input-> equal number of files required for fq1 and fq2"
#     input_str_fq1 = ",".join(fq1)
#     input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
#     input_str =  " ".join([input_str_fq1, input_str_fq2])

#     if fq1[0].endswith(".gz"):
#         readcmd = "--readFilesCommand zcat"
#     else:
#         readcmd = ""

#     outprefix = os.path.dirname(snakemake.output[0]) + "/"

#     shell(
#         "STAR "
#         "{extra} "
#         "--runThreadN {snakemake.threads} "
#         "--genomeDir {snakemake.params.index} "
#         "--readFilesIn {input_str} "
#         "{readcmd} "
#         "--outSAMtype BAM Unsorted "
#         "--outFileNamePrefix {outprefix} "
#         "--outStd Log "
#         "{log}")

# # https://snakemake-wrappers.readthedocs.io/en/0.17.4/wrappers/star/align.html
#     import os
#     from snakemake.shell import shell

#     extra = snakemake.params.get("extra", "")
#     log = snakemake.log_fmt_shell(stdout=True, stderr=True)

#     n = len(snakemake.input.sample)
#     assert n == 1 or n == 2, "input->sample must have 1 (single-end) or 2 (paired-end) elements."

#     if snakemake.input.sample[0].endswith(".gz"):
#         readcmd = "--readFilesCommand zcat"
#     else:
#         readcmd = ""


#     outprefix = os.path.dirname(snakemake.output[0]) + "/"


#     shell(
#         "STAR "
#         "{snakemake.params.extra} "
#         "--runThreadN {snakemake.threads} "
#         "--genomeDir {snakemake.params.index} "
#         "--readFilesIn {snakemake.input.sample} "
#         "{readcmd} "
#         "--outSAMtype BAM Unsorted "
#         "--outFileNamePrefix {outprefix} "
#         "--outStd Log "
#         "{log}")