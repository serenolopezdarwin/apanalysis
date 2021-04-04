import csv
import pickle as pkl
import sys


OVERLAP_PATH = "/net/shendure/vol1/home/sereno/projects/cell_clustering/nobackup/newannotations/data/overlapfiles/"


'''
This script is a subscript of matrixlengthandisoformanalysis.py
Its input and output are handled by matrixlengthandisoformanalysis.py
This script takes a single overlap file and calculates the distance of each read from the nearest PAS within 400 nt.
It returns a list of counts of each recorded distance, indexed from -400 to 400 relative to the PAS.
'''


def main():
    # input_path will be in the form of an overlap file the script will operate on.
    # pas_dataset will be in the form of a PAS dataset folder which contains files of annotated PAS loci.
    input_path = sys.argv[1]
    pas_dataset = sys.argv[2]
    reference_path = sys.argv[3]

    # This is a file of transcripts that have high 3'UTR-overlapping read counts, so we only calculate deviation on
    # highly-expressed genes.
    transcript_file_path = reference_path + "highly_expressed_transcripts.txt"
    # This folder is where all of the outputs from each of these subscripts will be stored for the head script to use.
    output_folder = reference_path + pas_dataset + "/pasmatrix/"

    # Loads the reference data used for the script. This includes information on genes and their PAS loci.
    with open(reference_path + pas_dataset + "/pas_data_dict.pkl", 'rb') as pas_file:
        pas_data_dict = pkl.load(pas_file)
    with open(reference_path + "names_by_id.pkl", 'rb') as names_in:
        cell_names = pkl.load(names_in)[0]

    input_file = input_path.replace(OVERLAP_PATH, reference_path + pas_dataset + "/assignedpas/").replace('.gz', '')
    # Makes a dictionary of the cells covered by this dataset and attaches the appropriate indices for our data.
    cell_chunk = {}
    with open(input_file, 'rt') as file_input:
        reader = csv.reader(file_input, delimiter='\t')
        for row in reader:
            cell = row[3]
            if cell in cell_chunk:
                continue
            elif cell not in cell_names:
                continue
            else:
                cell_chunk[cell] = cell_names.index(cell)

    # Reads each line of our input file and calculates expression level of each PAS of each gene.
    file_input = open(input_file, 'rt')
    reader = csv.reader(file_input, delimiter='\t')
    # This object will be built with a sub-dictionary for each gene with tuple keys for each cell, pas index.
    gene_dict = {}
    for row in reader:
        chrom = row[0]
        strand = row[5]
        gene = row[9]
        if chrom not in pas_data_dict or strand not in pas_data_dict[chrom] or gene not in pas_data_dict[chrom][strand]:
            continue

        cell = row[3]
        if cell not in cell_chunk:
            continue
        pas_locus = int(row[12])
        pas_length_idx = -1
        for pas_idx, pas_data in enumerate(pas_data_dict[chrom][strand][gene]):
            if pas_locus == pas_data[0]:
                pas_length_idx = pas_idx

        # If the pas length doesn't correspond to an existing PAS, skips this entry.
        if pas_length_idx == -1:
            continue

        # Checks if any of our data matches entered expression data, and adds it accordingly.
        if gene not in gene_dict:
            total_pas = [pas_data[1] for pas_data in pas_data_dict[chrom][strand][gene]]
            gene_dict[gene] = {'pas': total_pas}
        cell_idx = cell_chunk[cell]
        # Ignore the inspection, this works.
        if (cell_idx, pas_length_idx) not in gene_dict[gene]:
            gene_dict[gene][(cell_idx, pas_length_idx)] = 1
        else:
            gene_dict[gene][(cell_idx, pas_length_idx)] += 1

    file_input.close()

    # Stores our data as a dictionary to be merged into the global PAS data dictionary by the head script.
    # Our data is stored in the pasmatrix/ folder of our PAS dataset.
    out_file_path = output_folder + input_path.replace(".bed.gz", ".pkl").replace(OVERLAP_PATH, "")
    with open(out_file_path, 'wb') as pas_file_out:
        pkl.dump(gene_dict, pas_file_out)

    print("Output processed.")


if __name__ == "__main__":
    main()

