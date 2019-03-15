# Viper HIF workflow 

**Viper** is a [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html), aimed at performing the RNA-seq workflow for the project "HIF" in a reproducible, automated, and partially contained manner. It is implemented such that alternative or similar analysis can be added or removed. 

Viper consists of a `Snakefile` (`workflow/HIF_version_1.0/snakefile`), [`conda`](https://conda.io/docs/) environment files (`envs/*.yaml`), a configuration file (`workflow/HIF_version_1.0/config.yaml`), a set of `R` functions (`R/*R`), and a set of `R` scripts (`scripts/*.R`), to perform quality control, preprocessing, differential expression analysis, and functional annotation of RNA-seq data.

By default, the pipeline performs all the steps shown in the [diagram](img/report_2019_03_012_salmonAlignment_visualization.png) below. However advanced user, you can easly modiy `Snakefile` and `config.yaml` and/or "custom rules" to enble adtional functions. Currently, quantification with `Salmon` at the read-level or featureCounts can be enbled.

## Workflow graph
![DAG](img/report_2019_03_012_salmonAlignment_visualization.png)  

## Folder and File Structure 
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

## setup the VIPER workflow

Assuming that snakemake and conda are installed (and your system has the necessary libraries to compile R packages), you can use the following commands on a test dataset:

clone repo
```
git clonehttps://github.com/GrosseLab/ViperWF.git
```

make needed directories and copy files from `workflow/HIF_version_1.0`
```
mkdir data
mkdir data/qpcr
mkdir references
mkdir logs

cp ./viper/workflow/HIF_version_1.0/Snakefile
cp ./viper/workflow/HIF_version_1.0/config.yaml
cp ./viper/workflow/HIF_version_1.0/units.tsv
cp ./viper/workflow/HIF_version_1.0/samples.tsv
cp ./viper/workflow/HIF_version_1.0/cluster.json

cp ./viper/workflow/HIF_version_1.0/copy.csv ./data/qPCR/

cd ARMOR && snakemake --use-conda
```

# TODO: decription ### https://github.com/aerobatic/markdown-content/blob/master/docs/directory-structure.md  


0. Step - Step up folder  

git clone viper

