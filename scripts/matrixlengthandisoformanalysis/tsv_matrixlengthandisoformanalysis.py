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

    with open(REFERENCE_PATH + PAS_DATASET + "/cell_gene_expressions.pkl", 'rb') as exp_in:
        expressions = pkl.load(exp_in)
        cell_exps = expressions[0]
        gene_exps = expressions[1]
    # These are named after logarithmic cutoffs to be more distinct, i.e. 1 = 10^1 = 10, 4 = 10^4 = 10,000.
    approved_cells_1 = []
    approved_cells_2 = []
    approved_genes_3 = []
    approved_genes_4 = []
    for cell_id, cell_exp in enumerate(cell_exps):
        if cell_exp >= 100:
            approved_cells_1.append(cell_id)
            approved_cells_2.append(cell_id)
        elif cell_exp >= 10:
            approved_cells_1.append(cell_id)
    for gene_id, gene_exp in enumerate(gene_exps):
        if gene_exp >= 10000:
            approved_genes_3.append(gene_id)
            approved_genes_4.append(gene_id)
        elif gene_exp >= 1000:
            approved_genes_3.append(gene_id)

    h5_in_path = INPUT_FILE_PATH.replace(".bed.gz", ".h5") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/centered/")
    with h5.File(h5_in_path, 'r') as h5_in:
        cell_ids = list(h5_in['cells'])
        gene_ids = list(h5_in['genes'])
        utr_lengths = list(h5_in['utrs'])
        utr_bools = list(h5_in['length_bools'])
        cluster_lengths = list(h5_in['cluster_utrs'])
        trajectory_lengths = list(h5_in['traj_utrs'])
        subtrajectory_lengths = list(h5_in['subtraj_utrs'])
        age_lengths = list(h5_in['age_utrs'])

    with open(REFERENCE_PATH + "names_by_id.pkl", 'rb') as names_in:
        cell_names = pkl.load(names_in)[0]
    tsv_out_path = INPUT_FILE_PATH.replace(".bed.gz", ".tsv") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/tsv/")
    with open(tsv_out_path, 'wt') as cell_data_out:
        cell_utrs = []
        cell_bools = []
        cell_utrs_4 = []
        cell_utrs_3 = []
        cell_utrs_cluster = []
        cell_utrs_trajectory = []
        cell_utrs_subtrajectory = []
        cell_utrs_age = []
        for idx, cell_id in enumerate(cell_ids):
            gene_entry = gene_ids[idx]
            cell_utr = utr_lengths[idx]
            cell_bool = utr_bools[idx]
            cell_utr_cluster = cluster_lengths[idx]
            cell_utr_trajectory = trajectory_lengths[idx]
            cell_utr_subtrajectory = subtrajectory_lengths[idx]
            cell_utr_age = age_lengths[idx]
            cell_utrs.append(cell_utr)
            cell_bools.append(cell_bool)
            cell_utrs_cluster.append(cell_utr_cluster)
            cell_utrs_trajectory.append(cell_utr_trajectory)
            cell_utrs_subtrajectory.append(cell_utr_subtrajectory)
            cell_utrs_age.append(cell_utr_age)
            # Checks if the gene entry adheres to certain cutoffs, and adds it to separate 'high expression' utr lists.
            if gene_entry in approved_genes_4:
                cell_utrs_4.append(cell_utr)
                cell_utrs_3.append(cell_utr)
            elif gene_entry in approved_genes_3:
                cell_utrs_3.append(cell_utr)
            # Executes on the last cell group of the entire list or when a new cell group is on the next line.
            if idx + 1 == len(cell_ids) or cell_ids[idx + 1] != cell_id:
                cell_utr_mean = str(np.mean(cell_utrs))
                cell_ratio = str(np.mean(cell_bools))
                # Sets approved gene UTR means to 'NA' if cell has no reads from approved genes.
                # Otherwise this will set the approved gene UTR to the mean of only approved gene statistics.
                cell_utr_mean_4 = str(np.mean(cell_utrs_4)) if cell_utrs_4 else "NA"
                cell_utr_mean_3 = str(np.mean(cell_utrs_3)) if cell_utrs_3 else "NA"
                # Sets approved cell UTR means to 'NA' if the cells aren't in approved groups.
                if cell_id in approved_cells_2:
                    cell_utr_mean_2 = cell_utr_mean
                    cell_utr_mean_1 = cell_utr_mean
                elif cell_id in approved_cells_1:
                    cell_utr_mean_2 = "NA"
                    cell_utr_mean_1 = cell_utr_mean
                else:
                    cell_utr_mean_2 = "NA"
                    cell_utr_mean_1 = "NA"
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
                                  cell_utr_age_mean, cell_utr_mean_1, cell_utr_mean_2, cell_utr_mean_3, cell_utr_mean_4,
                                  cell_ratio, cell_data[20]]
                cell_data_str = '\t'.join(cell_data_used) + '\n'
                cell_data_out.write(cell_data_str)

                cell_utrs = []
                cell_bools = []
                cell_utrs_cluster = []
                cell_utrs_trajectory = []
                cell_utrs_subtrajectory = []
                cell_utrs_age = []
                cell_utrs_4 = []
                cell_utrs_3 = []

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

