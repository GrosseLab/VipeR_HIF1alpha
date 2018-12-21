

rule sambamba_sort:
    input:
        "{path}.bam"
    output:
        "{path}.sorted.bam"
    params:
        ""  # optional parameters
    threads: 20
    wrapper:
        "file:viper/wrapper/sambamba_v0.30.0/sort"        
        # "0.30.0/bio/sambamba/sort"

rule sambamba_index:
    input:
        "{path}.bam"
    output:
        "{path}.bam.bai"
    params:
        ""  # optional parameters
    threads: 20
    wrapper:
        "file:viper/wrapper/sambamba_Index"        

# rule sambamba_sort_index:
#     input:
#         "{path}.sorted.bam"
#     output:
#         "{path}.bam.bai"
#     params:
#         ""  # optional parameters
#     threads: 8
#     wrapper:
#         "file:viper/wrapper/sambamba_Index"         