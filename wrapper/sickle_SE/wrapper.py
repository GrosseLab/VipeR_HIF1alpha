__author__ = "Weinholdt Claus"
__copyright__ = "Copyright 2018, Weinholdt Claus"
__email__ = "claus.weinholdt@informatik.uni-halle.de."
__license__ = "MIT"

import os
from snakemake.shell import shell
from os import path

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("sickle se --version")

r1 = snakemake.input.get("r1")

outdir = path.dirname(snakemake.output.get('r1'))
shell("echo {r1}")
shell("pigz -dc -p {snakemake.threads} {r1} > {snakemake.output.r1tmp}")

shell(
	"(sickle se -f {snakemake.output.r1tmp} "
	"-o {snakemake.output.r1} "
	"-t {snakemake.params.qual_type} "
	"{extra}) {log}"
)

# __author__ = "Wibowo Arindrarto"
# __copyright__ = "Copyright 2016, Wibowo Arindrarto"
# __email__ = "bow@bow.web.id"
# __license__ = "BSD"

# from snakemake.shell import shell

# # Placeholder for optional parameters
# extra = snakemake.params.get("extra", "")
# log = snakemake.log_fmt_shell()

# shell(
#     "(sickle se -f {snakemake.input[0]} -o {snakemake.output[0]} "
#     "-t {snakemake.params.qual_type} {extra}) {log}"
#)
