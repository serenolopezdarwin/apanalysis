import csv
import gzip
import numpy as np
import pickle as pkl
import sys


INPUT_FILE_PATH = ""
OVERLAP_PATH = "/net/shendure/vol1/home/sereno/projects/cell_clustering/nobackup/newannotations/data/overlapfiles/"
PAS_DATASET = ""
CELL_DATA_DICT = {}
GENE_DATA_DICT = {}
REFERENCE_PATH = ""


'''
This script is a subscript of matrixrawandisoformanalysis.py
This script is set to be run in parallel. 
This script takes input in the form of a gzipped bed data file containing overlap data in a folder named "overlapfiles/"
This script will provide output in the form of three vectors corresponding to cell id, gene id, and utr length.
This script's output is based on raw, noncentered length calculations and therefore isn't very useful for analysis.
Instead, this script's output can be used to calculate averages across the entire dataset.
To test this script use:
>python raw_matrixlengthandisoformanalysis.py \
/net/shendure/vol1/home/sereno/projects/cell_clustering/nobackup/cellbed/overlapfiles/utr_overlaps_001.bed.gz raw

Test with:
python /net/shendure/vol1/home/sereno/projects/scripts/matrixlengthandisoformanalysis/\
raw_matrixlengthandisoformanalysis.py data/overlapfiles/utr_overlaps_###.bed.gz PAS_DATASET \
/net/shendure/vol10/projects/cell_clustering/nobackup/newannotations/
'''


def isoform_analysis():
    """Calculates PAS usage and exon coverage for each 3'UTR-overlapping read in our input data."""

    with open(REFERENCE_PATH + PAS_DATASET + "/pas_data_dict.pkl", 'rb') as pas_data_in:
        pas_data_dict = pkl.load(pas_data_in)
    # with open(REFERENCE_PATH + PAS_DATASET + "/pas_function.pkl", 'rb') as pas_function_in:
    #    pas_function = pkl.load(pas_function_in)
    pas_overlap_file_path = INPUT_FILE_PATH.replace(".bed.gz", ".bed")\
        .replace(OVERLAP_PATH, PAS_DATASET + "/assignedpas/")
    pas_overlap_file = open(pas_overlap_file_path, 'wt')

    outlier_count = 0
    file_input = gzip.open(INPUT_FILE_PATH, 'rt')
    reader = csv.reader(file_input, delimiter='\t')
    cell = 'first_cell'
    cell_age = 0
    cell_cluster = 0
    cell_trajectory = 0
    cell_subtrajectory = 0
    cell_dict = {}
    gene_dict = {}
    coverage_dict = {}
    locus_dict = {}
    for chrom in pas_data_dict:
        if chrom not in locus_dict:
            locus_dict[chrom] = {}
    for row in reader:
        chrom = row[0]
        strand = row[5]
        gene = row[9]
        if gene not in GENE_DATA_DICT or chrom not in pas_data_dict or strand not in pas_data_dict[chrom]\
                or gene not in pas_data_dict[chrom][strand]:
            continue

        new_cell = row[3]
        if new_cell not in CELL_DATA_DICT:
            continue
        # This makes our first cell id and cell dictionary.
        if cell == 'first cell':
            cell_age = CELL_DATA_DICT[new_cell][0]
            cell_cluster = CELL_DATA_DICT[new_cell][2]
            cell_trajectory = CELL_DATA_DICT[new_cell][8]
            cell_subtrajectory = CELL_DATA_DICT[new_cell][16]
        # This executes each time the script encounters a new cell.
        elif cell != new_cell:
            # First, we format the previous cell's data for output.
            for gene_id in cell_dict:
                cell_gene_median = np.median(cell_dict[gene_id])
                if gene_id not in gene_dict:
                    gene_dict[gene_id] = [[], {}, {}]
                if cell_cluster not in gene_dict[gene_id][1]:
                    gene_dict[gene_id][1][cell_cluster] = []
                if cell_trajectory not in gene_dict[gene_id][1]:
                    # Splitting this into two different indices of the head dict greatly simplifies later algorithms.
                    gene_dict[gene_id][1][cell_trajectory] = []
                if cell_subtrajectory not in gene_dict[gene_id][1]:
                    gene_dict[gene_id][1][cell_subtrajectory] = []
                # This looks like it should not be necessary. It is.
                if cell_trajectory not in gene_dict[gene_id][2]:
                    # This will hold our data for trajectories at certain ages.
                    gene_dict[gene_id][2][cell_trajectory] = {}
                if cell_age not in gene_dict[gene_id][2][cell_trajectory]:
                    gene_dict[gene_id][2][cell_trajectory][cell_age] = []
                if cell_age not in gene_dict[gene_id][1]:
                    gene_dict[gene_id][1][cell_age] = []
                gene_dict[gene_id][0].append(cell_gene_median)
                gene_dict[gene_id][1][cell_cluster].append(cell_gene_median)
                gene_dict[gene_id][1][cell_trajectory].append(cell_gene_median)
                gene_dict[gene_id][1][cell_age].append(cell_gene_median)
                gene_dict[gene_id][1][cell_subtrajectory].append(cell_gene_median)
                gene_dict[gene_id][2][cell_trajectory][cell_age].append(cell_gene_median)
            # Then, we reset our cell ID and cell dictionary for the next cell's entries.
            cell_age = CELL_DATA_DICT[new_cell][0]
            cell_cluster = CELL_DATA_DICT[new_cell][2]
            cell_trajectory = CELL_DATA_DICT[new_cell][8]
            cell_subtrajectory = CELL_DATA_DICT[new_cell][16]
            cell_dict = {}
        cell = new_cell

        # By default, we set the read as an outlier, which will change if we can attach it to an annotated isoform.
        # otherwise we scan forwards/backwards strand-wise
        gene_id = GENE_DATA_DICT[gene][0]
        pas_length = 30001
        pas_hits = []
        pas_counts = []
        pas_loci = []
        if strand == "+":
            locus = int(row[2])
            for pas_data in pas_data_dict[chrom][strand][gene]:
                pas = pas_data[0]
                pas_count = pas_data[3]
                pas_dist = locus - pas
                if -300 <= pas_dist <= 20:
                    pas_length = pas_data[1]
                    pas_hits.append((pas_dist, pas_length))
                    pas_loci.append(pas)
                    pas_counts.append(pas_count)
                else:
                    continue
        else:
            locus = int(row[1])
            for pas_data in pas_data_dict[chrom][strand][gene]:
                pas = pas_data[0]
                pas_count = pas_data[3]
                pas_dist = pas - locus
                if -300 <= pas_dist <= 20:
                    pas_length = pas_data[1]
                    pas_hits.append((pas_dist, pas_length))
                    pas_loci.append(pas)
                    pas_counts.append(pas_count)
                else:
                    continue

        # If no PAS overlaps are found, or the read's total length is an outlier, it is discarded.
        if pas_length > 30000:
            outlier_count = outlier_count + 1
            # Records the read's raw length for our coverage map.
            for exon in GENE_DATA_DICT[gene][1]:
                if exon[0] <= locus <= exon[1]:
                    if strand == "+":
                        raw_length = locus - exon[0] + exon[2]
                    else:
                        raw_length = exon[1] - locus + exon[2]
                    if gene not in coverage_dict:
                        coverage_dict[gene] = {'unfiltered': {raw_length: 1}}
                    elif raw_length in coverage_dict[gene]['unfiltered']:
                        coverage_dict[gene]['unfiltered'][raw_length] += 1
                    else:
                        coverage_dict[gene]['unfiltered'][raw_length] = 1
            continue
        else:
            # Writes the read and gene data to the pas assignment file.
            pas_overlap_file.write('\t'.join(row) + '\t')

        # If the read overlaps multiple PAS, we choose the nearest PAS and attach the read to it.
        if len(pas_hits) > 1:
            max_pas_idx = pas_counts.index(max(pas_counts))
            final_dist = pas_hits[max_pas_idx][0]
            final_length = pas_hits[max_pas_idx][1]
            locus = pas_loci[max_pas_idx]
            # Returns a list of tuples of pas hit indices and their corresponding probabilities
        #    pas_weights = [(idx, pas_function[((pas_hit[0] + 25) // 5)]) for idx, pas_hit in enumerate(pas_hits)]
            # Sorts the list of tuples by probabilities, takes the index of the entry with the greatest probability,
            # then references that index in pas_hits to retrieve the corresponding PAS length as our final length.
        #    final_pas = pas_hits[max(pas_weights, key=lambda t: t[1])[0]]
        #    final_dist = final_pas[0]
        #    final_length = final_pas[1]
        else:
            final_dist = pas_hits[0][0]
            final_length = pas_hits[0][1]
            locus = pas_loci[0]

        # Writes the PAS data to the pas assignment file.
        pas_overlap_file.write('\t'.join([str(locus), str(final_dist), str(final_length)]) + '\n')

        # Adds one to this locus' read assignment count
        if locus not in locus_dict[chrom]:
            locus_dict[chrom][locus] = 1
        else:
            locus_dict[chrom][locus] += 1

        raw_length = final_length + final_dist

        # Builds a dictionary of PAS lengths and their corresponding coverages, post-filtration.
        if gene not in coverage_dict:
            coverage_dict[gene] = {'filtered': {final_length: 1}, 'unfiltered': {raw_length: 1}}
        if 'filtered' not in coverage_dict[gene]:
            coverage_dict[gene] = {'filtered': {final_length: 1}}
        if 'unfiltered' not in coverage_dict[gene]:
            coverage_dict[gene] = {'unfiltered': {raw_length: 1}}
        if gene in coverage_dict and 'filtered' in coverage_dict[gene] and 'unfiltered' in coverage_dict[gene]:
            if final_length in coverage_dict[gene]['filtered']:
                coverage_dict[gene]['filtered'][final_length] += 1
            else:
                coverage_dict[gene]['filtered'][final_length] = 1
            if raw_length in coverage_dict[gene]['unfiltered']:
                coverage_dict[gene]['unfiltered'][raw_length] += 1
            else:
                coverage_dict[gene]['unfiltered'][raw_length] = 1

        if gene_id not in cell_dict:
            cell_dict[gene_id] = []
        cell_dict[gene_id].append(final_length)
    pas_overlap_file.close()

    # Executes on the last cell dataset of the file.
    for gene_id in cell_dict:
        cell_gene_median = np.median(cell_dict[gene_id])
        if gene_id not in gene_dict:
            gene_dict[gene_id] = [[], {}, {}]
        if cell_cluster not in gene_dict[gene_id][1]:
            gene_dict[gene_id][1][cell_cluster] = []
        if cell_trajectory not in gene_dict[gene_id][1]:
            gene_dict[gene_id][1][cell_trajectory] = []
        if cell_subtrajectory not in gene_dict[gene_id][1]:
            gene_dict[gene_id][1][cell_subtrajectory] = []
        if cell_trajectory not in gene_dict[gene_id][2]:
            gene_dict[gene_id][2][cell_trajectory] = {}
        if cell_age not in gene_dict[gene_id][2][cell_trajectory]:
            gene_dict[gene_id][2][cell_trajectory][cell_age] = []
        if cell_age not in gene_dict[gene_id][1]:
            gene_dict[gene_id][1][cell_age] = []
        gene_dict[gene_id][0].append(cell_gene_median)
        gene_dict[gene_id][1][cell_cluster].append(cell_gene_median)
        gene_dict[gene_id][1][cell_trajectory].append(cell_gene_median)
        gene_dict[gene_id][1][cell_age].append(cell_gene_median)
        gene_dict[gene_id][2][cell_trajectory][cell_age].append(cell_gene_median)
        gene_dict[gene_id][1][cell_subtrajectory].append(cell_gene_median)

    file_input.close()

    print(outlier_count)

    return gene_dict, coverage_dict, locus_dict


"""def calculate_threepseq_coverage(coverage_dict_build):
    #Formats and adds the 3pseq data to our coverage dictionary.

    threepseq_in_file_path = INPUT_FILE_PATH.replace("utr_overlaps", "3pseq")\
        .replace(OVERLAP_PATH, "data/3pseq/")
    file_input = gzip.open(threepseq_in_file_path, 'rt')
    reader = csv.reader(file_input, delimiter='\t')

    for row in reader:
        gene = row[9]
        if gene not in GENE_DATA_DICT:
            continue

        strand = row[5]
        if strand == "+":
            locus = int(row[2])
        else:
            locus = int(row[1])

        coverage = float(row[3])
        for exon in GENE_DATA_DICT[gene][1]:
            if exon[0] <= locus <= exon[1]:
                if strand == "+":
                    length = locus - exon[0] + exon[2]
                else:
                    length = exon[1] - locus + exon[2]
                if gene not in coverage_dict_build:
                    coverage_dict_build[gene] = {'3pseq': {length: coverage}}
                elif '3pseq' not in coverage_dict_build[gene]:
                    coverage_dict_build[gene]['3pseq'] = {length: coverage}
                elif length in coverage_dict_build[gene]['3pseq']:
                    coverage_dict_build[gene]['3pseq'][length] += coverage
                else:
                    coverage_dict_build[gene]['3pseq'][length] = coverage

    file_input.close()

    coverage_dict = coverage_dict_build

    return coverage_dict"""


def main():
    global INPUT_FILE_PATH
    INPUT_FILE_PATH = sys.argv[1]
    global PAS_DATASET
    PAS_DATASET = sys.argv[2]
    global REFERENCE_PATH
    REFERENCE_PATH = sys.argv[3]

    with open(REFERENCE_PATH + "cell_data_dict.pkl", 'rb') as cell_data_in:
        global CELL_DATA_DICT
        CELL_DATA_DICT = pkl.load(cell_data_in)
    with open(REFERENCE_PATH + "gene_data_dict.pkl", 'rb') as gene_data_in:
        global GENE_DATA_DICT
        GENE_DATA_DICT = pkl.load(gene_data_in)

    gene_dict, coverage_dict_build, locus_dict = isoform_analysis()

    #coverage_dict = calculate_threepseq_coverage(coverage_dict_build)

    # Our output is stored in a pkl file in the raw/ folder of our PAS dataset folder.
    print(INPUT_FILE_PATH)
    print(OVERLAP_PATH)
    print(REFERENCE_PATH)
    raw_out_file_path = INPUT_FILE_PATH.replace(".bed.gz", ".pkl")\
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/raw/")
    with open(raw_out_file_path, 'wb') as raw_out_file:
        pkl.dump(gene_dict, raw_out_file)
    print(raw_out_file_path)

    #coverage_out_file_path = INPUT_FILE_PATH.replace(".bed.gz", ".pkl")\
    #    .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/coverage/")
    #with open(coverage_out_file_path, 'wb') as coverage_out_file:
    #    pkl.dump(coverage_dict, coverage_out_file)

    locus_file_out_path = INPUT_FILE_PATH.replace(".bed.gz", ".pkl")\
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/pasloci/")
    with open(locus_file_out_path, 'wb') as locus_out_file:
        pkl.dump(locus_dict, locus_out_file)
    print(locus_file_out_path)

if __name__ == "__main__":
    main()

