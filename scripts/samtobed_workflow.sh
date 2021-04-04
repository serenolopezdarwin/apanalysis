#!/bin/bash

printf -v jobid "%03d" ${SGE_TASK_ID}
overlappath="data/overlapfiles/utr_overlaps_$jobid.bed"

# Converts our sam read data to more useable bed format data. 
bedfile="data/cellbed/data_chunk_$jobid.bed"
while read f; do
    fileid=$(basename $f .sam)
    sam2bed <$f | awk -v id=$fileid 'BEGIN{OFS=FS="\t"}{print $1, $2, $3, id, "1", $6}'
done < temp$jobid | bedtools sort > $bedfile

# Creates files for the overlaps between our new bed file and our gene annotations.
bedtools intersect -s -a $bedfile -b full_overlap.bed -s -wa -wb | cut -f 10 > temp$jobid.transcripts
bedtools intersect -s -a $bedfile -b 3utr_overlap.bed -wa -wb > $overlappath
sort -k4,4 -o $overlappath $overlappath
gzip $overlappath

# Processes the chunk's corresponding 3pseq bulk data file.
threepfile="data/3pseq/3pseq_$jobid.bed"
awk 'BEGIN{OFS=FS="\t"}{print $1, $2, $3, $5, "1", $4}' $threepfile | bedtools sort | \
bedtools intersect -s -a stdin -b 3utr_overlap.bed -wa -wb > temp2_$jobid
mv temp2_$jobid $threepfile
gzip $threepfile

echo "Output processed."
