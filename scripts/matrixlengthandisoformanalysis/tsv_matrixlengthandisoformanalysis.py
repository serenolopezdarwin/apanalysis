import h5py as h5
import numpy as np
import pickle as pkl
import sys


INPUT_FILE_PATH = ""
OVERLAP_PATH = "data/overlapfiles/"
PAS_DATASET = ""
CELL_DATA_DICT = {}
REFERENCE_PATH = ""


'''
This script is a subscript of matrixrawandisoformanalysis.py
This script is set to be run in parallel. 
This script takes input in the form of a hdf5 file containing a sparse cell-by-gene utr length matrix in a folder named
 "centered/" in the PAS dataset directory.
This script will provide output in the form of tsv files .
This script's output is formatted for easy plotting.
'''


def generate_cell_tsv():
    """This function is essentially the entire script. It reads, our raw cell by gene data and calculates various
    dataset means based on filtering of data points by certain pre-determined cell or gene expression levels."""

    h5_in_path = INPUT_FILE_PATH.replace(".bed.gz", ".h5") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/centered/")
    with h5.File(h5_in_path, 'r') as h5_in:
        cell_ids = list(h5_in['cells'])
        utr_lengths = list(h5_in['utrs'])
        cluster_lengths = list(h5_in['cluster_utrs'])
        trajectory_lengths = list(h5_in['traj_utrs'])
        subtrajectory_lengths = list(h5_in['subtraj_utrs'])
        age_lengths = list(h5_in['age_utrs'])

    with open(REFERENCE_PATH + "names_by_id.pkl", 'rb') as names_in:
        cell_names = pkl.load(names_in)[0]
    tsv_out_path = INPUT_FILE_PATH.replace(".bed.gz", ".tsv") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/tsv/")
    with open(tsv_out_path, 'wt') as cell_data_out:
        cell_count = 0
        cell_utrs = []
        cell_utrs_cluster = []
        cell_utrs_trajectory = []
        cell_utrs_subtrajectory = []
        cell_utrs_age = []
        for idx, cell_id in enumerate(cell_ids):
            cell_count += 1
            cell_utr = utr_lengths[idx]
            cell_utr_cluster = cluster_lengths[idx]
            cell_utr_trajectory = trajectory_lengths[idx]
            cell_utr_subtrajectory = subtrajectory_lengths[idx]
            cell_utr_age = age_lengths[idx]
            cell_utrs.append(cell_utr)
            cell_utrs_cluster.append(cell_utr_cluster)
            cell_utrs_trajectory.append(cell_utr_trajectory)
            cell_utrs_subtrajectory.append(cell_utr_subtrajectory)
            cell_utrs_age.append(cell_utr_age)
            # Executes on the last cell group of the entire list or when a new cell group is on the next line.
            if idx + 1 == len(cell_ids) or cell_ids[idx + 1] != cell_id:
                cell_utr_mean = str(np.mean(cell_utrs))
                # Sets approved gene UTR means to 'NA' if cell has no reads from approved genes.
                # Otherwise this will set the approved gene UTR to the mean of only approved gene statistics.
                # Sets approved cell UTR means to 'NA' if the cells aren't in approved groups.
                cell_utr_cluster_mean = str(np.mean(cell_utrs_cluster))
                cell_utr_trajectory_mean = str(np.mean(cell_utrs_trajectory))
                cell_utr_subtrajectory_mean = str(np.mean(cell_utrs_subtrajectory))
                cell_utr_age_mean = str(np.mean(cell_utrs_age))
                cell_name = cell_names[cell_id]
                cell_data = CELL_DATA_DICT[cell_name]
                cell_age = cell_data[0]
                cell_subcluster = cell_data[2] + "." + cell_data[5]
                cell_data_used = [cell_data[2], cell_data[3], cell_data[4], cell_subcluster, cell_data[6], cell_data[7],
                                  cell_data[8], cell_data[9], cell_data[10], cell_data[11], cell_data[16],
                                  cell_data[13], cell_data[14], cell_data[15], cell_age, cell_utr_mean,
                                  cell_utr_cluster_mean, cell_utr_trajectory_mean, cell_utr_subtrajectory_mean,
                                  cell_utr_age_mean, cell_data[20], cell_count]
                cell_data_str = '\t'.join(cell_data_used) + '\n'
                cell_data_out.write(cell_data_str)
                # Resets cell data for next line.
                cell_utrs = []
                cell_utrs_cluster = []
                cell_utrs_trajectory = []
                cell_utrs_subtrajectory = []
                cell_utrs_age = []
                cell_count = 0

    print("Cell tsv generated!")


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

    generate_cell_tsv()


if __name__ == "__main__":
    main()
