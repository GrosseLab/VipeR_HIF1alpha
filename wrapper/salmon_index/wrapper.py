__author__ = "Weinholdt Claus"
__copyright__ = "Copyright 2018, Weinholdt Claus"
__email__ = "claus.weinholdt@informatik.uni-halle.de."
__license__ = "MIT"

import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("salmon index -t {snakemake.input} -i {snakemake.output} "
      " --threads {snakemake.threads} {extra} {log}" )

# """Snakemake wrapper for Salmon Index."""

# __author__ = "Tessa Pierce"
# __copyright__ = "Copyright 2018, Tessa Pierce"
# __email__ = "ntpierce@gmail.com"
# __license__ = "MIT"

# from snakemake.shell import shell

# log = snakemake.log_fmt_shell(stdout=True, stderr=True)
# extra = snakemake.params.get("extra", "")

# shell("salmon index -t {snakemake.input} -i {snakemake.output} "
#       " --threads {snakemake.threads} {extra} {log}" )