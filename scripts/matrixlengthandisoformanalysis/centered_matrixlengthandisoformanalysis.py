import csv
import gzip
import h5py as h5
import numpy as np
import pickle as pkl
import sys


INPUT_FILE_PATH = ""
OVERLAP_PATH = "data/overlapfiles/"
PAS_DATASET = ""
CELL_DATA_DICT = {}
GENE_DATA_DICT = {}
GENE_AVERAGE_DICT = {}
REFERENCE_PATH = ""


'''
This script is a subscript of matrixrawandisoformanalysis.py
This script is set to be run in parallel. 
This script takes input in the form of a gzipped bed data file containing overlap data in a folder named "overlapfiles/"
This script will provide output in the form of four vectors corresponding to cell id, gene id, utr length, and counts.
This script will also provide cell and gene expression level counts for this chunk of data.
This script's output is formatted for final data analysis.
'''


def isoform_analysis():
    with open(REFERENCE_PATH + PAS_DATASET + "/pas_data_dict.pkl", 'rb') as pas_data_in:
        pas_data_dict = pkl.load(pas_data_in)
    # with open(REFERENCE_PATH + PAS_DATASET + "/pas_function.pkl", 'rb') as pas_function_in:
    #     pas_function = pkl.load(pas_function_in)
    with open(REFERENCE_PATH + "/names_by_id.pkl", 'rb') as names_in:
        names = pkl.load(names_in)
        pas_names = names[3]
        total_genes = len(names[1])
        total_pas = len(pas_names)

    file_input = gzip.open(INPUT_FILE_PATH, 'rt')
    reader = csv.reader(file_input, delimiter='\t')
    cell_ids = []
    gene_ids = []
    utr_lengths = []
    length_bools = []
    cluster_lengths = []
    trajectory_lengths = []
    subtrajectory_lengths = []
    age_lengths = []
    counts = []
    cell_counts = []
    cell_expression_level = 0
    gene_counts = [0] * total_genes
    cell = 'first_cell'
    cell_id = 0
    cell_cluster = 0
    cell_trajectory = 0
    cell_subtrajectory = 0
    cell_age = 0
    missing_averages = 0
    cell_dict = {}

    ages = ['9.5', '10.5', '11.5', '12.5', '13.5']
    clusters = ['Cardiac_muscle_lineages', 'Cholinergic_neurons', 'Chondroctye_progenitors',
                'Chondrocytes_and_osteoblasts', 'Connective_tissue_progenitors', 'Definitive_erythroid_lineage',
                'Early_mesenchyme', 'Endothelial_cells', 'Ependymal_cell', 'Epithelial_cells', 'Excitatory_neurons',
                'Granule_neurons', 'Hepatocytes', 'Inhibitory_interneurons', 'Inhibitory_neuron_progenitors',
                'Inhibitory_neurons', 'Intermediate_Mesoderm', 'Isthmic_organizer_cells', 'Jaw_and_tooth_progenitors',
                'Lens', 'Limb_mesenchyme', 'Megakaryocytes', 'Melanocytes', 'Myocytes', 'Neural_progenitor_cells',
                'Neural_Tube', 'Neutrophils', 'Notochord_cells', 'Oligodendrocyte_Progenitors', 'Osteoblasts',
                'Postmitotic_premature_neurons', 'Premature_oligodendrocyte', 'Primitive_erythroid_lineage',
                'Radial_glia', 'Schwann_cell_precursor', 'Sensory_neurons', 'Stromal_cells', 'White_blood_cells']
    cell_age_idx = 0
    cell_cluster_idx = 0
    gene_array = np.zeros((total_genes, 5))
    count_array = np.zeros((total_genes, 5))
    cluster_gene_array = np.zeros((total_genes, 38))
    cluster_count_array = np.zeros((total_genes, 38))
    trajectory_array_dict = {}

    age_pas_array = np.zeros((total_pas, 5))
    cluster_pas_array = np.zeros((total_pas, 38))

    for row in reader:
        chrom = row[0]
        strand = row[5]
        gene = row[9]
        if gene not in GENE_DATA_DICT or chrom not in pas_data_dict or strand not in pas_data_dict[chrom] \
                or gene not in pas_data_dict[chrom][strand]:
            continue

        new_cell = row[3]
        if new_cell not in CELL_DATA_DICT:
            continue
        # This makes our first cell id and cell dictionary.
        if cell == 'first cell':
            cell_id = CELL_DATA_DICT[new_cell][19]
            cell_age = CELL_DATA_DICT[new_cell][0]
            cell_age_idx = ages.index(cell_age)
            cell_cluster = CELL_DATA_DICT[new_cell][2]
            if cell_cluster != "NA":
                cell_cluster_idx = clusters.index(cell_cluster)
            cell_trajectory = CELL_DATA_DICT[new_cell][8]
            cell_subtrajectory = CELL_DATA_DICT[new_cell][16]
            cell_expression_level = 0
        # This executes each time the script encounters a new cell.
        elif cell != new_cell:
            # First, we format the previous cell's data for output.
            for gene_id in cell_dict:
                # Marks expression of this gene in a single cell for later filtration
                gene_counts[gene_id] += 1
                # Processing utr data by centering it to different dataset means.
                cell_gene_median = np.median(cell_dict[gene_id])
                try:
                    gene_average = GENE_AVERAGE_DICT[gene_id][0]
                    cell_gene_median_centered = cell_gene_median - gene_average
                    cluster_average = GENE_AVERAGE_DICT[gene_id][1][cell_cluster]
                    cell_gene_median_cluster = cell_gene_median - cluster_average
                    trajectory_average = GENE_AVERAGE_DICT[gene_id][1][cell_trajectory]
                    cell_gene_median_trajectory = cell_gene_median - trajectory_average
                    subtrajectory_average = GENE_AVERAGE_DICT[gene_id][1][cell_subtrajectory]
                    cell_gene_median_subtrajectory = cell_gene_median - subtrajectory_average
                    age_average = GENE_AVERAGE_DICT[gene_id][1][cell_age]
                except KeyError:
                    missing_averages += 1
                    continue
                cell_gene_median_age = cell_gene_median - age_average
                cell_ids.append(cell_id)
                gene_ids.append(gene_id)
                utr_lengths.append(cell_gene_median_centered)
                if cell_gene_median_centered >= 0:
                    length_bools.append(1)
                else:
                    length_bools.append(0)
                cluster_lengths.append(cell_gene_median_cluster)
                trajectory_lengths.append(cell_gene_median_trajectory)
                subtrajectory_lengths.append(cell_gene_median_subtrajectory)
                age_lengths.append(cell_gene_median_age)
                counts.append(len(cell_dict[gene_id]))
                cell_gene_count = len(cell_dict[gene_id])
                counts.append(cell_gene_count)
                cell_expression_level += cell_gene_count
                # Here we build our gene data arrays. We do some very basic mean re-adjustment if there is already data
                # for this gene at this age (from another cell).
                old_count = count_array[gene_id, cell_age_idx]
                old_mean = gene_array[gene_id, cell_age_idx]
                # We only add one to our data count for each cell to normalize the data to gene expression levels.
                total_count = old_count + 1
                count_array[gene_id, cell_age_idx] = total_count
                new_mean = (old_mean * old_count + cell_gene_median_centered) / total_count
                gene_array[gene_id, cell_age_idx] = new_mean
                if cell_cluster != "NA":
                    cluster_old_count = cluster_count_array[gene_id, cell_cluster_idx]
                    cluster_old_mean = cluster_gene_array[gene_id, cell_cluster_idx]
                    cluster_total_count = old_count + 1
                    cluster_count_array[gene_id, cell_cluster_idx] = cluster_total_count
                    cluster_new_mean = (cluster_old_mean * cluster_old_count + cell_gene_median_centered) / total_count
                    cluster_gene_array[gene_id, cell_cluster_idx] = cluster_new_mean
                if cell_trajectory not in trajectory_array_dict:
                    # Initiates a count array and a gene array the first time we encounter a certain trajectory.
                    # Index 0 is the gene array, index 1 is the count array.
                    trajectory_array_dict[cell_trajectory] = [np.zeros((total_genes, 5)), np.zeros((total_genes, 5))]
                trajectory_gene_array = trajectory_array_dict[cell_trajectory][0]
                trajectory_count_array = trajectory_array_dict[cell_trajectory][1]
                old_count = trajectory_count_array[gene_id, cell_age_idx]
                old_mean = trajectory_gene_array[gene_id, cell_age_idx]
                total_count = old_count + 1
                trajectory_array_dict[cell_trajectory][1][gene_id, cell_age_idx] = total_count
                # We use the trajectory-adjusted median to calculate our mean values
                new_mean = (old_mean * old_count + cell_gene_median_trajectory) / total_count
                trajectory_array_dict[cell_trajectory][0][gene_id, cell_age_idx] = new_mean

            # Then, we reset our cell data for the next cell's entries.
            cell_counts.append((cell_expression_level, cell_id))
            cell_id = CELL_DATA_DICT[new_cell][19]
            cell_age = CELL_DATA_DICT[new_cell][0]
            cell_age_idx = ages.index(cell_age)
            cell_cluster = CELL_DATA_DICT[new_cell][2]
            if cell_cluster != "NA":
                cell_cluster_idx = clusters.index(cell_cluster)
            cell_trajectory = CELL_DATA_DICT[new_cell][8]
            cell_subtrajectory = CELL_DATA_DICT[new_cell][16]
            cell_dict = {}
            cell_expression_level = 0
        cell = new_cell

        # By default, we set the read as an outlier, which will change if we can attach it to an annotated isoform.
        gene_id = GENE_DATA_DICT[gene][0]
        pas_length = 30001
        pas_hits = []
        pas_counts = []
        if strand == "+":
            locus = int(row[2])
            for pas_idx, pas_data in enumerate(pas_data_dict[chrom][strand][gene]):
                pas = pas_data[0]
                pas_count = pas_data[3]
                pas_dist = locus - pas
                if -300 <= pas_dist <= 20:
                    pas_length = pas_data[1]
                    pas_hits.append((pas_length, pas_idx))
                    pas_counts.append(pas_count)
                else:
                    continue
        else:
            locus = int(row[1])
            for pas_idx, pas_data in enumerate(pas_data_dict[chrom][strand][gene]):
                pas = pas_data[0]
                pas_count = pas_data[3]
                pas_dist = pas - locus
                if -300 <= pas_dist <= 20:
                    pas_length = pas_data[1]
                    pas_hits.append((pas_length, pas_idx))
                    pas_counts.append(pas_count)
                else:
                    continue

        # If no PAS overlaps are found, or the read's total length is an outlier, it is discarded.
        if pas_length > 30000:
            continue

        # If the read overlaps multiple PAS, we choose the nearest PAS and attach the read to it.
        if len(pas_hits) > 1:
            final_hit_idx = pas_counts.index(max(pas_counts))
            final_length = pas_hits[final_hit_idx][0]
            final_pas_idx = pas_hits[final_hit_idx][1]
            # Returns a list of tuples of pas hit indices and their corresponding probabilities based on our pas
            # function and the read's distance from each pas it overlaps.
            # pas_weights = [(idx, pas_function[((pas_hit[0] + 25) // 5)]) for idx, pas_hit in enumerate(pas_hits)]
            # Sorts the list of tuples by probabilities, takes the index of the entry with the greatest probability,
            # then references that index in pas_hits to retrieve the corresponding PAS length and idx.
            # final_hit_idx = max(pas_weights, key=lambda t: t[1])[0]
            # final_length = pas_hits[final_hit_idx][1]
            # final_pas_idx = pas_hits[final_hit_idx][2]
        else:
            final_length = pas_hits[0][0]
            final_pas_idx = pas_hits[0][1]

        if gene_id not in cell_dict:
            cell_dict[gene_id] = []
        cell_dict[gene_id].append(final_length)

        total_pas_idx = pas_names.index(gene + "." + str(final_pas_idx))
        age_pas_array[(total_pas_idx, cell_age_idx)] += 1
        cluster_pas_array[(total_pas_idx, cell_cluster_idx)] += 1

    # This executes on the last cell dataset of the file.
    for gene_id in cell_dict:
        gene_counts[gene_id] += 1
        cell_gene_median = np.median(cell_dict[gene_id])
        try:
            gene_average = GENE_AVERAGE_DICT[gene_id][0]
            cell_gene_median_centered = cell_gene_median - gene_average
            cluster_average = GENE_AVERAGE_DICT[gene_id][1][cell_cluster]
            cell_gene_median_cluster = cell_gene_median - cluster_average
            trajectory_average = GENE_AVERAGE_DICT[gene_id][1][cell_trajectory]
            cell_gene_median_trajectory = cell_gene_median - trajectory_average
            subtrajectory_average = GENE_AVERAGE_DICT[gene_id][1][cell_subtrajectory]
            cell_gene_median_subtrajectory = cell_gene_median - subtrajectory_average
            age_average = GENE_AVERAGE_DICT[gene_id][1][cell_age]
        except KeyError:
            missing_averages += 1
            continue
        if cell_gene_median_centered >= 0:
            length_bools.append(1)
        else:
            length_bools.append(0)
        cell_gene_median_age = cell_gene_median - age_average
        cell_ids.append(cell_id)
        gene_ids.append(gene_id)
        utr_lengths.append(cell_gene_median_centered)
        cluster_lengths.append(cell_gene_median_cluster)
        trajectory_lengths.append(cell_gene_median_trajectory)
        subtrajectory_lengths.append(cell_gene_median_subtrajectory)
        age_lengths.append(cell_gene_median_age)
        cell_gene_count = len(cell_dict[gene_id])
        counts.append(cell_gene_count)
        cell_expression_level += cell_gene_count
        # Here we build our gene data arrays. We do some very basic mean re-adjustment if there is already data
        # for this gene at this age (from another cell).
        old_count = count_array[gene_id, cell_age_idx]
        old_mean = gene_array[gene_id, cell_age_idx]
        # We only add one to our data count for each cell to normalize the data to gene expression levels.
        total_count = old_count + 1
        count_array[gene_id, cell_age_idx] = total_count
        new_mean = (old_mean * old_count + cell_gene_median_centered) / total_count
        gene_array[gene_id, cell_age_idx] = new_mean
        if cell_cluster != "NA":
            cluster_old_count = cluster_count_array[gene_id, cell_cluster_idx]
            cluster_old_mean = cluster_gene_array[gene_id, cell_cluster_idx]
            cluster_total_count = old_count + 1
            cluster_count_array[gene_id, cell_cluster_idx] = cluster_total_count
            cluster_new_mean = (cluster_old_mean * cluster_old_count + cell_gene_median_centered) / total_count
            cluster_gene_array[gene_id, cell_cluster_idx] = cluster_new_mean
        if cell_trajectory not in trajectory_array_dict:
            # Initiates a count array and a gene array the first time we encounter a certain trajectory.
            # Index 0 is the gene array, index 1 is the count array.
            trajectory_array_dict[cell_trajectory] = [np.zeros((total_genes, 5)), np.zeros((total_genes, 5))]
        trajectory_gene_array = trajectory_array_dict[cell_trajectory][0]
        trajectory_count_array = trajectory_array_dict[cell_trajectory][1]
        old_count = trajectory_count_array[gene_id, cell_age_idx]
        old_mean = trajectory_gene_array[gene_id, cell_age_idx]
        total_count = old_count + 1
        trajectory_array_dict[cell_trajectory][1][gene_id, cell_age_idx] = total_count
        # We use the trajectory-adjusted median to calculate our mean values
        new_mean = (old_mean * old_count + cell_gene_median_trajectory) / total_count
        trajectory_array_dict[cell_trajectory][0][gene_id, cell_age_idx] = new_mean
    cell_counts.append((cell_expression_level, cell_id))

    return cell_ids, gene_ids, utr_lengths, cluster_lengths, trajectory_lengths, age_lengths, counts, \
        subtrajectory_lengths, gene_array, count_array, trajectory_array_dict, cell_counts, gene_counts, \
        cluster_gene_array, cluster_count_array, age_pas_array, cluster_pas_array, length_bools, missing_averages


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
    with open(REFERENCE_PATH + PAS_DATASET + "/gene_average_dict.pkl", 'rb') as gene_average_in:
        global GENE_AVERAGE_DICT
        GENE_AVERAGE_DICT = pkl.load(gene_average_in)

    cell_ids, gene_ids, utr_lengths, cluster_lengths, trajectory_lengths, age_lengths, counts, \
     subtraj_lengths, gene_array, count_array, trajectory_array_dict, cell_counts, gene_counts, \
     cluster_gene_array, cluster_count_array, age_pas_array, cluster_pas_array, length_bools, \
     missing_averages = isoform_analysis()

    array_out_path = INPUT_FILE_PATH.replace(".bed.gz", ".pkl") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/heatmap/")
    with open(array_out_path, 'wb') as array_out_file:
        pkl.dump([gene_array, count_array, trajectory_array_dict,
                  cluster_gene_array, cluster_count_array], array_out_file)

    pas_stats_out_path = INPUT_FILE_PATH.replace(".bed.gz", ".pkl") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/pasarrays/")
    with open(pas_stats_out_path, 'wb') as pas_stats_out_file:
        pkl.dump([age_pas_array, cluster_pas_array], pas_stats_out_file)

    # Our output is stored in a h5 file in the length/ folder of our PAS dataset folder.
    h5_out_path = INPUT_FILE_PATH.replace(".bed.gz", ".h5") \
        .replace(OVERLAP_PATH, REFERENCE_PATH + PAS_DATASET + "/centered/")
    with h5.File(h5_out_path, "w") as out_file:
        out_file.create_dataset("cells", data=cell_ids)
        out_file.create_dataset("genes", data=gene_ids)
        out_file.create_dataset("utrs", data=utr_lengths)
        out_file.create_dataset("cluster_utrs", data=cluster_lengths)
        out_file.create_dataset("traj_utrs", data=trajectory_lengths)
        out_file.create_dataset("age_utrs", data=age_lengths)
        out_file.create_dataset("subtraj_utrs", data=subtraj_lengths)
        out_file.create_dataset("counts", data=counts)
        out_file.create_dataset("cell_counts", data=cell_counts)
        out_file.create_dataset("gene_counts", data=gene_counts)
        out_file.create_dataset("length_bools", data=length_bools)

    print(missing_averages)


if __name__ == "__main__":
    main()

