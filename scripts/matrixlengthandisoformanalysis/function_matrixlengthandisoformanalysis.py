import csv
import gzip
import numpy as np
import pickle as pkl
import sys


OVERLAP_PATH = "data/overlapfiles/"


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

    transcript_file_path = reference_path + "detected_genes_parsedoverlaps_1transcript_per_gene_OKIDs.txt"
    # This folder is where all of the outputs from each of these subscripts will be stored for the head script to use.
    output_folder = reference_path + pas_dataset + "/function/"

    # Loads the reference data used for the script. This includes information on genes and their PAS loci.
    with open(transcript_file_path, 'rt') as transcript_file:
        approved_transcripts = [line.rstrip('\n') for line in transcript_file]
    with open(reference_path + pas_dataset + "/pas_data_dict.pkl", 'rb') as pas_file:
        pas_data_dict = pkl.load(pas_file)

    # Reads each line of our input file and calculates which PAS the read belongs to and how far it is from that PAS.
    input_file = sys.argv[1]
    # Use 161 as the length for -400 to 400
    hit_list = np.array([0] * 86)
    file_input = gzip.open(input_file, 'rt')
    reader = csv.reader(file_input, delimiter='\t')
    for row in reader:
        chrom = row[0]
        strand = row[5]
        gene = row[9]
        if gene not in approved_transcripts or chrom not in pas_data_dict or strand not in pas_data_dict[chrom] \
                or gene not in pas_data_dict[chrom][strand]:
            continue

        # This checks if our read is within 400 NT of one or more PAS and calculates how far it is from each PAS.
        pas_hits = []
        pas_counts = []
        if strand == "+":
            locus = int(row[2])
            for pas_data in pas_data_dict[chrom][strand][gene]:
                pas = pas_data[0]
                pas_count = pas_data[3]
                pas_dist = locus - pas
                if -400 <= pas_dist <= 25:
                    pas_hits.append(pas_dist)
                    pas_counts.append(pas_count)
                else:
                    continue
        else:
            locus = int(row[1])
            for pas_data in reversed(pas_data_dict[chrom][strand][gene]):
                pas = pas_data[0]
                pas_count = pas_data[3]
                pas_dist = pas - locus
                if -400 <= pas_dist <= 25:
                    pas_hits.append(pas_dist)
                    pas_counts.append(pas_count)
                else:
                    continue

        # In the case of a read being within range of multiple PAS, we eliminate it so as not to create noise.
        # These reads will be modified using the function this analysis creates later on, so we do not want them
        # contaminating the original function.
        if len(pas_hits) == 0:
            continue
        elif len(pas_hits) == 1:
            hit_dist = pas_hits[0]
        # If the read overlaps multiple PAS, takes the one with the greatest read count.
        else:
            pas_max_idx = pas_counts.index(max(pas_counts))
            hit_dist = pas_hits[pas_max_idx]
        hit_index = (400 + hit_dist) // 5
        hit_list[hit_index] = hit_list[hit_index] + 1

    # Stores our data as a numpy array to be concatenated with other numpy arrays by the head script.
    # Our data is stored in the function/ folder of our PAS dataset.
    hit_array = np.array(hit_list)
    out_file_path = output_folder + input_path.replace(".bed.gz", ".pkl").replace(OVERLAP_PATH, "")
    with open(out_file_path, 'wb') as pas_file_out:
        pkl.dump(hit_array, pas_file_out)

    print("Output processed.")


if __name__ == "__main__":
    main()

