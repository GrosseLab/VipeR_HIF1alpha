.. This workflow performs differential expression analysis on paired-end RNA-seq data.
.. After adapter removal with `Cutadapt <http://cutadapt.readthedocs.io>`_ and quality filtering with `sickle <https://github.com/najoshi/sickle>`_, reads were mapped with `STAR <https://github.com/alexdobin/STAR>`_ to the humane genome (GRCh38.82). 
.. The gene counts were generated with  `featureCounts <http://subread.sourceforge.net>`_ in three different modes -- 'unique counting', 'fraction counting', and 'multiple mapped counting'.
.. Further, reads were mapped with `salmon <https://github.com/COMBINE-lab/salmon>`_ and STAR to the humane transcriptome (GRCh38.82) and transcript counts were generated with salmon. These transcript counts were summarized to gene counts with `tximport <https://github.com/mikelove/tximport>`_. 
.. Integrated normalization and differential expression analysis was conducted with `edegR <https://bioconductor.org/packages/release/bioc/html/edgeR.html>`_.

This workflow performs differential expression analysis on paired-end RNA-seq data.
After adapter removal with `Cutadapt <http://cutadapt.readthedocs.io>`_ and quality filtering with `sickle <https://github.com/najoshi/sickle>`_, reads were mapped with `STAR <https://github.com/alexdobin/STAR>`_ to the humane genome (GRCh38.82) and transcript counts were quantified with `salmon <https://github.com/COMBINE-lab/salmon>`_ . 
These transcript counts were summarized to gene counts with `tximport <https://github.com/mikelove/tximport>`_. 
Integrated normalization and differential expression analysis was conducted with `edegR <https://bioconductor.org/packages/release/bioc/html/edgeR.html>`_. 
Further  we used the Database for Annotation, Visualization and Integrated Discovery (`DAVID v6_8  <https://david.ncifcrf.gov/content.jsp?file=citation.html>`_) for a functional annotation of the differential expressed genes.
