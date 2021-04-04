import csv
import pickle as pkl
import sys


INPUT_FILE_PATH = ""
OVERLAP_PATH = "data/overlapfiles/"
PAS_DATASET = ""
REFERENCE_PATH = ""


'''
This script is a subscript of matrixrawandisoformanalysis.py
This script is set to be run in parallel. 
This script takes input in the form of a hdf5 file containing a sparse cell-by-gene utr length matrix in a folder named
 "centered/" in the PAS dataset directory.
This script will provide output in the form of tsv files .
This script's output is formatted for easy plotting.
'''


def count_transcripts():
    """This function is essentially the entire script. It reads, our raw cell by gene data and calculates various
    dataset means based on filtering of data points by certain pre-determined cell or gene expression levels."""

    transcript_counts = {}
    pas_overlap_in_path = INPUT_FILE_PATH.replace(".gz", "") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/assignedpas/")
    with open(pas_overlap_in_path, 'rt') as overlaps_in:
        reader = csv.reader(overlaps_in, delimiter='\t')
        for row in reader:
            transcript = row[9]
            if transcript not in transcript_counts:
                transcript_counts[transcript] = 0
            else:
                transcript_counts[transcript] += 1

    transcript_counts_out_path = INPUT_FILE_PATH.replace(".bed.gz", ".pkl") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/count/")
    with open(transcript_counts_out_path, 'wb') as transcript_counts_out:
        pkl.dump(transcript_counts, transcript_counts_out)

    print("Transcripts counted!")


def main():
    global INPUT_FILE_PATH
    INPUT_FILE_PATH = sys.argv[1]
    global PAS_DATASET
    PAS_DATASET = sys.argv[2]
    global REFERENCE_PATH
    REFERENCE_PATH = sys.argv[3]

    count_transcripts()


if __name__ == "__main__":
    main()

