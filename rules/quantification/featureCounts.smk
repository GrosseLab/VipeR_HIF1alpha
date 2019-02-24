# rule feature_counts:
#     input:
#         # Required input. Recommend using wildcards for sample names,
#         # e.g. {sample,SRR[0-9]+}
#         bam = "results/mapping/star/{ref}/{sample}-{unit}/Aligned.out.sorted.bam",             # ['{sample}.bam'],
#         annotation =lambda wildcards: config["ref"][wildcards.ref]["annotation"]                     #'annotation/gencode.v28.gtf.gz'
#     output:
#         # Required output.
#         'results/quantification/featureCounts/{ref}/{sample}-{unit}.feature_counts.result'
#     params:
#         # Optional parameters. Omit if unneeded.
#         extra = '',
#         countparam = '',
#         # Format of provided annotation file. SAF or GTF.
#         annotation_format = 'GTF',
#         # GTF-specific options.
#         feature_type = lambda wildcards: config["ref"][wildcards.ref]["sjdbGTFfeatureExon"], ##'exon',  # Feature type in GTF file.
#         attribute_type = 'gene_id',  # Attribute type in GTF file.
#         # Minimum number of overlapping bases in a read that is required for read assignment.
#         # [default: 1]
#         min_overlap = 1,
#     threads: 24
#     wrapper:
#         "file:viper/wrapper/featureCounts"
#         #'http://dohlee-bio.info:9193/feature-counts'

rule featureCounts:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        bam = "results/mapping/star/{ref}/{sample}-{unit}/Aligned.out.sorted.bam",             # ['{sample}.bam'],
        annotation =lambda wildcards: config["ref"][wildcards.ref]["annotation"]                     #'annotation/gencode.v28.gtf.gz'
    output:
        # Required output.
        res='results/quantification/featureCounts/{ref}/{readtype}/{ctype}/{sample}-{unit}.feature_counts.result',
        summary='results/quantification/featureCounts/{ref}/{readtype}/{ctype}/{sample}-{unit}.feature_counts.result.summary'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        countparam = lambda wildcards: config["featureCounts"][wildcards.readtype][wildcards.ctype],
        # Format of provided annotation file. SAF or GTF.
        annotation_format = 'GTF',
        # GTF-specific options.
        feature_type = lambda wildcards: config["ref"][wildcards.ref]["sjdbGTFfeatureExon"], ##'exon',  # Feature type in GTF file.
        attribute_type = 'gene_id',  # Attribute type in GTF file.
        # Minimum number of overlapping bases in a read that is required for read assignment.
        # [default: 1]
        min_overlap = 1,
    threads: 2
    wrapper:
        "file:viper/wrapper/featureCounts"
        #'http://dohlee-bio.info:9193/feature-counts'