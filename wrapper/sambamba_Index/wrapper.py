__author__ = "Weinholdt Claus"
__copyright__ = "Copyright 2018, Weinholdt Claus"
__email__ = "claus.weinholdt@informatik.uni-halle.de."
__license__ = "MIT"


from snakemake.shell import shell

shell("sambamba index {snakemake.params} -t {snakemake.threads} {snakemake.input[0]} {snakemake.output[0]}")
