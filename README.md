# Analyzing APA With sci-RNA-seq Data 

This repository is intended to accompany our publication, primarily to enhance the reproducibility of our results. For more information please refer to:

Agarwal V, Lopez-Darwin S., Kelley D., Shendure J. [The landscape of alternative polyadenylation in single cells of the developing mouse embryo.](https://www.biorxiv.org/content/10.1101/2021.01.21.427498v1)


This pipeline can be used on sci-rna-seq datasets in BAM format in order to quantify and plot 3'UTR usage at up to single-cell and single-gene resolution.

If you find our code or predictions to be helpful for your work, please cite the paper above.


# Dependencies for running entire pipeline:
* Python3 modules: glob2, h5py, matplotlib, numpy, pandas, pylab, scipy, seaborn, subprocess, matplotlib_venn, statsmodels

* R libraries: ggplot2, ggseqlogo, gplots, LSD, pheatmap, plotly, plyr, dplyr, processx, RColorBrewer, reshape2, viridis

* [The MEME Suite](http://meme-suite.org/doc/download.html?man_type=web)

* [UCSC tools](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads) installation

* [BEDTools](https://github.com/arq5x/bedtools2/releases)

# Instructions for use


Users are advised to read the code closely and modify commented pieces as appropriate to acquire
desired output for your environment. For example, you will need to download all of the additional
R library and Python module dependencies for the code to work. This being said, if you find crucial
files are missing, making the code unusable, or if you identify a major problem in the code, please
raise a Github issue.

Directory management is generally self-contained but note that some directories, such as the script path and initial data directory path must be defined in the head script workflow.sh in order to properly function.
