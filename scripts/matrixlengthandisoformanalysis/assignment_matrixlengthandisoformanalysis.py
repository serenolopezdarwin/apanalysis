import csv
import gzip
import pickle as pkl
import sys


OVERLAP_PATH = "/net/shendure/vol1/home/sereno/projects/cell_clustering/nobackup/newannotations/data/overlapfiles"


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
    output_folder = reference_path + pas_dataset + "/assignment/"

    # Loads the reference data used for the script. This includes information on genes and their PAS loci.
    with open(transcript_file_path, 'rt') as transcript_file:
        approved_transcripts = [line.rstrip('\n') for line in transcript_file]
    with open(reference_path + pas_dataset + "/pas_data_dict.pkl", 'rb') as pas_file:
        pas_data_dict = pkl.load(pas_file)

    # Reads each line of our input file and calculates which PAS the read belongs to and how far it is from that PAS.
    # A dictionary identical to pas_data_dict but contains a read count for each PAS.
    pas_counts = pas_data_dict
    file_input = gzip.open(input_path, 'rt')
    reader = csv.reader(file_input, delimiter='\t')
    for row in reader:
        chrom = row[0]
        strand = row[5]
        gene = row[9]
        if gene not in approved_transcripts or chrom not in pas_data_dict or strand not in pas_data_dict[chrom] \
                or gene not in pas_data_dict[chrom][strand]:
            continue

        # This checks if our read is within 400 NT of one or more PAS and calculates how far it is from each PAS.
        if strand == "+":
            locus = int(row[2])
            for pas_idx, pas_data in enumerate(pas_data_dict[chrom][strand][gene]):
                pas = pas_data[0]
                pas_dist = locus - pas
                if -300 <= pas_dist <= 20:
                    # Adds one to the pas read count if the pas is in the allowed window.
                    pas_counts[chrom][strand][gene][pas_idx][3] += 1
                else:
                    continue
        else:
            locus = int(row[1])
            for pas_idx, pas_data in enumerate(pas_data_dict[chrom][strand][gene]):
                pas = pas_data[0]
                pas_dist = pas - locus
                if -300 <= pas_dist <= 20:
                    pas_counts[chrom][strand][gene][pas_idx][3] += 1
                else:
                    continue

    # Stores our data as a dictionary to be merged into the global PAS data dictionary by the head script.
    # Our data is stored in the assignment/ folder of our PAS dataset.
    out_file_path = output_folder + input_path.replace(".bed.gz", ".pkl").replace(OVERLAP_PATH, "")
    with open(out_file_path, 'wb') as pas_file_out:
        pkl.dump(pas_counts, pas_file_out)

    print("Output processed.")


if __name__ == "__main__":
    main()

