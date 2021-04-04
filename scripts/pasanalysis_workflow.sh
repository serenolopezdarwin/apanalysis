#!/bin/bash

# These two lines (4, 5) is replaced by the head script.
dataset="cutoff"
wd="/net/shendure/vol1/home/sereno/projects/cell_clustering/nobackup/newannotations"
scriptpath=/net/shendure/vol1/home/sereno/projects/scripts/

Rscript ${scriptpath}differentialpasanalysis.R $dataset $wd ${SGE_TASK_ID}
