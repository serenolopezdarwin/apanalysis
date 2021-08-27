#!/bin/bash

    
# This script requires one argument, being the name of the dataset to be analyzed (no '/' needed.)
# This analysis requires the following files: 
# 1. Several folders full of reads in SAM format 
# 2. A set of gene annotations in gtf format named "gencode.vM25.annotation.gtf"
# 3. A set of PAS annotations in bed format named "PAS.bed" with six fields corresponding to
## Chromosome, start, stop, name, score, strand
# 4. A job file for processing the raw data named "samtobed_workflow.sh"
# 5. A text file of baseline expression levels of mm10 genes named "mouse_median_expr.txt"
# 6. A python script to process and analyze our data, named 'matrixlengthandisoformanalysis.py'
## This script requires the following subscripts: 
## a. function_matrixlengthandisoformanalysis.py
## b. raw_matrixlengthandisoformanalysis.py
## c. count_matrixlengthandisoformanalysis.py
## d. centered_matrixlengthandisoformanalysis.py
## e. tsv_matrixlengthandisoformanalysis.py
## This script has the following library dependencies:
## glob2, h5py, matplotlib, matplotlib_venn, numpy, pandas, pylab, scipy, and seaborn.
# 7. A csv file of cell annotations named "cell_annotations.csv"
# 8. A fasta file containing the mm10 sequence named 'mm10.fa' and its alignment file named "mm10.fa.fai"
# 9. A file of 3pseq bulk data in bed format named "3pseq.bed"
# 10. A R script to plot some of our formatted data, named 'generatedataplots.R'
## This script has the following library dependencies:
## ggplot2, gplots, LSD, pheatmap, plotly, plyr, RColorBrewer, reshape2, and viridis.

# This analysis requires the cli tools:
# 1. bedtools
# 2. gtf2bed
# 3. sge
# 4. bedops


# If no dataset is entered, everything is directed to the 'default' dataset folder.
dataset=${1:-default}
scriptpath=/net/shendure/vol1/home/sereno/projects/scripts/
# Point this at your python path.
pythonpath=~sereno/software/anaconda2/bin/python

# Makes some folders we will need to hold our formatted data.
if [ ! -d data ]
then
    mkdir data
fi
if [ ! -d data/cellbed ]
then
    mkdir data/cellbed
fi
if [ ! -d data/overlapfiles ]
then
    mkdir data/overlapfiles
fi
if [ ! -d data/pasloci ]
then
    mkdir data/pasloci
fi
if [ ! -d data/3pseq ]
then
    mkdir data/3pseq
fi
if [ ! -d $dataset ]
then
    mkdir $dataset
fi
if [ ! -d $dataset/pasoverlaps ]
then
    mkdir $dataset/pasoverlaps
fi
if [ ! -d $dataset/pasloci ]
then
    mkdir $dataset/pasloci
fi
if [ ! -d $dataset/figures ]
then
    mkdir $dataset/figures
fi
if [ ! -d $dataset/figures/transcriptviolin/ ]
then
    mkdir $dataset/figures/transcriptviolin
fi
if [ ! -d $dataset/figures/transcriptviolin/ages ]
then
    mkdir $dataset/figures/transcriptviolin/ages
fi
if [ ! -d $dataset/figures/transcriptviolin/clusters ]
then
    mkdir $dataset/figures/transcriptviolin/clusters
fi
if [ ! -d $dataset/figures/transcriptviolin/trajectories ]
then
    mkdir $dataset/figures/transcriptviolin/trajectories
fi
if [ ! -d $dataset/assignedpas ]
then
    mkdir $dataset/assignedpas
fi
if [ ! -d $dataset/significantpas ]
then
    mkdir $dataset/significantpas
fi

# Formats our gene annotations
if [ ! -e "3utr_overlap.bed" ]
then
    # Takes entire gene body of each gene, with introns included.
    sort -k9,9 -k4,4n gene_annotations.gtf | \
    awk 'BEGIN{OFS=FS="\t"} \ 
    {if (!a[$9]++) {curr=$0;$0=prev; print $5, $9, "1", $7; $0=curr; print $1, $4}} \
    {prev=$0} END{print $5, $9, "1", $7}' | awk NR\>1 | paste -d"\t" - - | bedtools sort > full_overlap.bed
    gtf2bed < gene_annotations.gtf | grep 'three_prime_UTR' | \
    awk 'BEGIN{OFS=FS="\t"}{print $1, $2, $3, $9, "1", $5}' | bedtools sort > 3utr_overlap.bed
fi

# Submits parallel jobs which will process our raw sam files.
if [ ! -e "transcript_counts_overlapfiles.txt" ]
then
    echo "Processing raw data..."
    # Splits our files into 4000-file chunks.
    # Point this at your raw data.
    raw_data_path='/net/shendure/vol1/home/sereno/projects/cell_clustering/data'
    # If you have more than 4M files, I envy you, but you must use -a >3 instead.
    find $raw_data_path/ -name '*.sam' | split --numeric-suffixes=1 -a 3 -l 4000 - temp
    shopt -s nullglob
    numfiles=(temp*)
    numfiles=${#numfiles[@]}

    # Splits our 3pseq file into the same number of chunks as our other data, for later concurrent parallelism.
    chunks=$numfiles
    linetotal=$(wc -l < '3pseq.bed')
    ((chunklines = ($linetotal + $chunks - 1) / $chunks))
    split --numeric-suffixes=1 --additional-suffix=.bed -a 3 -l $chunklines 3pseq.bed data/3pseq/3pseq_ 

    # Submits a separate job for each of our data chunks.
    # Point this wherever the job file is.
    jobfile='/net/shendure/vol1/home/sereno/projects/scripts/samtobed_workflow.sh'
    # If you have less than 20 data chunks, you must use -tc <20 instead.
    qsub -t 1-$numfiles -tc 20 -N samtobed bash $jobfile

    # Concatenates our results and counts up the occurences of each transcript for later analysis.
    # It is faster to do this as a multi-cat than as a for loop, hence the use of temp files.
    cat *.transcripts | sort | uniq -c > transcript_counts_original.txt
    rm temp*
    for f in data/overlapfiles/*.bed.gz; do
        zcat $f | cut -f 10 > $f.temp
    done 
    # 
    cat *.bed.gz.temp | sort | uniq -c > transcript_counts_overlapfiles.txt
    rm *temp

    echo "Raw data processed."
else
    echo "Raw data already processed."
fi

# If the PAS annotation superset isn't provided, this will generate it from our files.
# To create your own PAS annotation set simply use the following bash command on your full annotations:
## bedtools sort -i {YOURPAS.bed} | \
## bedtools intersect -s -a stdin -b 3utr_overlap.bed -wa -wb > $dataset/PAS.bed
if [ -e "PAS.bed" ]
then
    echo "PAS annotations found."
    mv PAS.bed $dataset/PAS.bed
elif [ -e "$dataset/PAS.bed" ]
then
    echo "Pas annotations found."
else
    echo "Generating PAS annotations..."
    # Rename these 'original' files to whatever PAS input files you're using.
    # Fixes our umdnj annotations by extending them, fixing chromosome annotation, and overlapping with 
    # our 3' UTR annotations.
    awk 'BEGIN{OFS=FS="\t"}{ if ($6=="+") {$2=$3-10} else {$3=$2+10}}{print $1, $2, $3, "PolyADBv3", "1", $6}' \
    umdnj_PAS_original.bed | sed -e 's/^/chr/' | bedtools sort | \
    bedtools intersect -a stdin -b 3utr_overlap.bed -s > umdnj_extend.temp.bed

    #cat .bed umdnj_extend.temp.bed | \
    #perl -ne 'chomp; @a=split/\t/; if($a[5] eq "+"){ $a[1]=$a[2]-50; $a[2]+=50; } else{ $a[2]=$a[1]+50; $a[1]-=50; }; \
    #print join("\t", @a),"\n";' | \
    #bedtools getfasta -fi /net/shendure/vol10/nobackup/shared/genomes/mus_musculus/mm10.fa -bed - \
    #-tab -name -s -fo /dev/stdout | Rscript ggseqlogo.R

    # Extracts the most annotated PAS coordinate from the Unibas mm10 annotations, fixes chromosome annotation,
    # extends reads, and overlaps with vikram's extended UTR annotations.
    awk 'BEGIN{OFS=FS="\t"}{split($4,a,":"); if ($6=="+") {$3=a[2]; $2=$3-10} else {$2=a[2]; $3=$2+10}} \
    {print $1, $2, $3, "PolyASite2", "1", $6}' atlas.clusters.2.0.GRCm38.96.bed | sed -e 's/^/chr/' | \
    bedtools sort | bedtools intersect -s -a stdin -b 3utr_overlap.bed > unibas_extend.temp.bed

    # Extracts the end coordinate of the most downstream 3'UTR of gencode annotations and extends by 10 NT.
    awk '{if ($3=="UTR") print}' gencode.vM25.annotation.gtf | \
    awk 'BEGIN{OFS=FS="\t"}{split($9,a,";"); split(a[2],b,"\""); \
    $6=b[2]}{print $1, $4, $5, $9, "1", $7}' | \
    tac | awk -F"\t" '!a[$4]++' | \
    awk 'BEGIN{OFS=FS="\t"}{if ($6=="+") {$2=$3-10} else {$3=$2+10}}{print $1, $2, $3, "Gencode", $5, $6}' | \
    bedtools sort | bedtools intersect -s -a stdin -b 3utr_overlap.bed > gencode_extend.temp.bed

    # Repeatedly searches through existing PAS annotations to find unannotated ones from each dataset, 
    # appending the new annotations each time onto the existing ones. 
    # This creates a union of all three of our input datasets we created above, and ensures uniqueness.
    bedtools intersect -s -v -a umdnj_extend.temp.bed -b unibas_extend.temp.bed | \
    cat - unibas_extend.temp.bed | bedtools sort | \
    bedtools intersect -s -v -a gencode_extend.temp.bed -b stdin | cat - umdnj_unibas.temp.bed | \
    uniq | bedtools sort | uniq | bedtools intersect -s -a stdin -b 3utr_overlap.bed -wa -wb | \
    sed '/zero_length_insertion/d' > $dataset/PAS.bed

    # Generates unique sites from each of our PAS source datasets and builds some plots from them.
    bedtools intersect -s -v -a unibas_extend.temp.bed -b gencode_extend.temp.bed | \
    cat - gencode_extend.temp.bed | bedtools sort | \
    bedtools intersect -s -v -a umdnj_extend.temp.bed -b stdin | \
    awk 'BEGIN{OFS=FS="\t"}{if ($6=="+") {$2=$3-1} else {$3=$2+1}}{print $1, $2, $3, $4, $5, $6}' | \
    bedtools sort > umdnj_unique.temp.bed
    cat umdnj_unique.bed | perl -ne 'chomp; @a=split/\t/; if($a[5] eq "+"){ $a[1]=$a[2]-50; $a[2]+=50; }
    else{ $a[2]=$a[1]+50; $a[1]-=50; }; print join("\t", @a),"\n";' | \
    bedtools getfasta -fi /net/shendure/vol10/nobackup/shared/genomes/mus_musculus/mm10.fa \
    -bed - -tab -name -s -fo /dev/stdout | Rscript ggseqlogo.R
    mv figs/Fig1D.PFM.all.original.100K.pdf figs/polyadb_pfm.pdf

    bedtools intersect -s -v -a umdnj_extend.temp.bed -b gencode_extend.temp.bed | \
    cat - gencode_extend.temp.bed | bedtools sort | \
    bedtools intersect -s -v -a unibas_extend.temp.bed -b stdin | \
    awk 'BEGIN{OFS=FS="\t"}{if ($6=="+") {$2=$3-1} else {$3=$2+1}}{print $1, $2, $3, $4, $5, $6}' | \
    bedtools sort > unibas_unique.temp.bed
    cat unibas_unique.temp.bed | perl -ne 'chomp; @a=split/\t/; if($a[5] eq "+"){ $a[1]=$a[2]-50; $a[2]+=50; }
    else{ $a[2]=$a[1]+50; $a[1]-=50; }; print join("\t", @a),"\n";' | \
    bedtools getfasta -fi /net/shendure/vol10/nobackup/shared/genomes/mus_musculus/mm10.fa \
    -bed - -tab -name -s -fo /dev/stdout | Rscript ggseqlogo.R
    mv figs/Fig1D.PFM.all.original.100K.pdf figs/polyasite_pfm.pdf

    bedtools intersect -s -v -a umdnj_extend.temp.bed -b unibas_extend.temp.bed | \
    cat - unibas_extend.temp.bed | bedtools sort | \
    bedtools intersect -s -v -a gencode_extend.temp.bed -b stdin | \
    awk 'BEGIN{OFS=FS="\t"}{if ($6=="+") {$2=$3-1} else {$3=$2+1}}{print $1, $2, $3, $4, $5, $6}' | \
    bedtools sort > gencode_unique.temp.bed
    cat gencode_unique.temp.bed | perl -ne 'chomp; @a=split/\t/; if($a[5] eq "+"){ $a[1]=$a[2]-50; $a[2]+=50; }
    else{ $a[2]=$a[1]+50; $a[1]-=50; }; print join("\t", @a),"\n";' | \
    bedtools getfasta -fi /net/shendure/vol10/nobackup/shared/genomes/mus_musculus/mm10.fa \
    -bed - -tab -name -s -fo /dev/stdout | Rscript ggseqlogo.R
    mv figs/Fig1D.PFM.all.original.100K.pdf figs/gencode_pfm.pdf

    rm *.temp.*

    echo "PAS annotations generated."
fi
if [ -e $dataset/PAS_ext.bed ]; 
then
    echo "PAS counts processed."
else
    # Creates a new version of our PAS file that extends each PAS' range by 800 to find overlapping reads.
    paspath="$dataset/PAS.bed"
    awk 'BEGIN{OFS=FS="\t"}{if ($6=="+") {$2=$3-400;$3=$2+800} else {$3=$2+400;$2=$3-800}}
    {print $1, $2, $3, $4, $5, $6}' $paspath | uniq > $dataset/PAS_ext.bed

    shopt -s nullglob
    numfiles=(data/overlapfiles/*.bed.gz)
    numfiles=${#numfiles[@]}

    # Runs the pas-overlapping script in parallel on our read data
    jobfile='/net/shendure/vol1/home/sereno/projects/scripts/pasoverlap_workflow.sh'
    sed -i "6s|.*|dataset=\"$dataset\"|" $jobfile
    qsub -t 1-$numfiles -tc 20 -N pasoverlap bash $jobfile
fi


# This python script analyzes and formats our data. Its output is detailed in its docstring.
$pythonpath ${scriptpath}matrixlengthandisoformanalysis.py $dataset $PWD

# It is faster to do this here than in the python script
cat $dataset/tsv/*.tsv > $dataset/cell_data.tsv
# This will submit 100 parallel jobs to analyze differential PAS usage.
# Notice: This fragment has been deprecated in favor of marking differential PAS inside of the python script.
# jobfile='/net/shendure/vol1/home/sereno/projects/scripts/pasanalysis_workflow.sh'
# sed -i "4s|.*|dataset=\"$dataset\"|" $jobfile
# sed -i "5s|.*|wd=\"$PWD\"|" $jobfile
# qsub -t 1-100 -tc 10 -N pasanalysis bash $jobfile

# This R script generates data plots. Its output is detailed in the script.
Rscript ${scriptpath}generatedataplots.R $dataset $PWD
