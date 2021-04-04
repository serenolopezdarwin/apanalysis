#!/bin/bash

printf -v jobid "%03d" ${SGE_TASK_ID}
overlappath="data/overlapfiles/utr_overlaps_$jobid.bed.gz"
# This line (6) is replaced by the head script.
dataset="newtx"
paspath="$dataset/PAS_ext.bed"
pasoverlappath="$dataset/pasoverlaps/pas_overlaps_$jobid.bed"

# Creates files for the overlaps between our bed file and our PAS annotations.
zcat $overlappath | awk 'BEGIN{OFS=FS="\t"}{ if ($6=="+") {$2=$3} else {$3=$2}}{print $1, $2, $3, $4, $5, $6}' | \
bedtools sort | bedtools intersect -s -a stdin -b $paspath -c | \
awk 'BEGIN{OFS=FS="\t"}{ if ($7==1) {print $1, $2, $3, $4, $5, $6}}' | \
bedtools intersect -s -a stdin -b $paspath -wa -wb | \
awk 'BEGIN{OFS=FS="\t"}{if ($6=="+") {$13=($3-$8)-400} else {$13=($9-$2)-400}}
{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' > $pasoverlappath

echo "Output processed."
