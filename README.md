# ViperWF
Snakemake workflow "viper"

[Snakemake workflow "viper" on Slurm](./report_2018_12_21.html)


Folder and File Structure 
============================

Here is the basic suggested skeleton for your project repo:

```bash
  .
  ├── data
  │   ├── qPCR 	            # qPRCR raw data
  │   └ *.fastq.gz 	        # all 'fastq.gz'-files from !...!
  │
  ├── references
  │   └── hg38 	    				      # all data from Homo_sapiens.GRCh38.82
  │   ├ Homo_sapiens.GRCh38.82.gtf 	    				          # annotation
  │   ├ Homo_sapiens.GRCh38.dna.primary_assembly.fa 	        # genome sequence 
  │   └ Homo_sapiens.GRCh38.82.EXON.fa 	    				      # exon sequence of all transcript of GTF
  │	
  ├── logs
  │
  ├── viper 	    				      # Github repository 
  │   ├── report 	    				      # Snakemake report definition
  │   ├── wrapper 	    				      # Snakemake wrapper
  │   ├── rules 	    				      # Snakemake rules
  │   ├── scripts 	    				      # Snakemake scripts
  │   ├── workflow 	    				      # Snakemake final workflows
  │   │	  └ HIF_version_1.0 	    				      #
  │   ├── R 	    				      # R functions needed to tun analysis   
  │   └── man 	    				      # R functions manual
  │
  ├── Snakefile 	    				      # file from ./viper/workflow/HIF_version_1.0
  ├── config.yaml 	    				      # file from ./viper/workflow/HIF_version_1.0
  ├── units.tsv 	    				      # file from ./viper/workflow/HIF_version_1.0
  ├── samples.tsv 	    				      # file from ./viper/workflow/HIF_version_1.0
  └── cluster.json 	    				      # file from ./viper/workflow/HIF_version_1.0
```

# TODO: decription ### https://github.com/aerobatic/markdown-content/blob/master/docs/directory-structure.md  


0. Step - Step up folder  

git clone viper
mkdir data
mkdir references
mkdir logs

cp ./viper/workflow/HIF_version_1.0/Snakefile
cp ./viper/workflow/HIF_version_1.0/config.yaml
cp ./viper/workflow/HIF_version_1.0/units.tsv
cp ./viper/workflow/HIF_version_1.0/samples.tsv
cp ./viper/workflow/HIF_version_1.0/cluster.json