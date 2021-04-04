import csv
import glob2
import h5py as h5
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle as pkl
import pylab as pl
import random as rd
import re
import scipy.stats as sps
import scipy.cluster.hierarchy as sch
import seaborn as sns
import subprocess
import sys
import time
# noinspection PyUnresolvedReferences
from matplotlib_venn import venn3

OVERLAP_PATH = "data/overlapfiles/"
TRANSCRIPT_FILE_PATH = "detected_genes_parsedoverlaps_1transcript_per_gene_OKIDs.txt"
SCRIPT_PATH = "/net/shendure/vol1/home/sereno/projects/scripts/matrixlengthandisoformanalysis/"
PAS_DATASET = sys.argv[1] if len(sys.argv) > 1 else "raw"
APPROVED_TRANSCRIPTS = []
GENE_NAMES = []
AGES = ['9.5', '10.5', '11.5', '12.5', '13.5']
CELL_TYPES = ['Cardiac_muscle_lineages', 'Cholinergic_neurons', 'Chondroctye_progenitors',
              'Chondrocytes_and_osteoblasts', 'Connective_tissue_progenitors', 'Definitive_erythroid_lineage',
              'Early_mesenchyme', 'Endothelial_cells', 'Ependymal_cell', 'Epithelial_cells', 'Excitatory_neurons',
              'Granule_neurons', 'Hepatocytes', 'Inhibitory_interneurons', 'Inhibitory_neuron_progenitors',
              'Inhibitory_neurons', 'Intermediate_Mesoderm', 'Isthmic_organizer_cells', 'Jaw_and_tooth_progenitors',
              'Lens', 'Limb_mesenchyme', 'Megakaryocytes', 'Melanocytes', 'Myocytes', 'Neural_progenitor_cells',
              'Neural_Tube', 'Neutrophils', 'Notochord_cells', 'Oligodendrocyte_Progenitors', 'Osteoblasts',
              'Postmitotic_premature_neurons', 'Premature_oligodendrocyte', 'Primitive_erythroid_lineage',
              'Radial_glia', 'Schwann_cell_precursor', 'Sensory_neurons', 'Stromal_cells', 'White_blood_cells']

'''
This script calculates utr deviation values for each individual cell
This script's inputs are handled by the workflow manager script "workflow.sh"
This script gives the following outputs:
1. A pickled dictionary of cell annotation data named "cell_data_dict.pkl"
2. A pickled list of lists of cell, transcript, and gene names ordered by their ids named "names_by_id.pkl"
3. A pickled dictionary of gene read data named "gene_data_dict.pkl"
4. A pickled dictionary of gene PAS annotation data named "pas_data_dict.pkl" in our dataset folder.
5. A pickled list of weights corresponding to NT distances from annotated PAS sites, named "pas_function.pkl"
6. A pickled dictionary of gene utr length averages, named "gene_average_dict.pkl" in our dataset folder.
7. An array of gene data grouped by cell ages named "gene_age_array.pkl" in our dataset folder.
8. A TSV file of read counts for each gene at each round of filtration in our dataset folder.
9. An hdf5 file containing vectors corresponding to sparse matrix vectors. These vectors are as follows:
   cell id, gene id, utr length (centered to global gene mean), and data count for that cell x gene coordinate.
10. Several folders corresponding to our generated data, separated into our original data chunks. Each data chunk
    corresponds to 4000 cells. This data includes cell ids, gene ids, and utr lengths centered to a variety of factors 
    such as dataset mean utr for the gene in certain clusters, trajectories, or ages.
11. A variety of plots. They are contained in the datasets's figures/ folder. These plots are as follows:
    a. A venn diagram of the shared annotation of PAS between our source datasets named "pas_dataset_overlap.pdf"
    b. A lineplot of the distribution of our reads relative to our PAS dataset names "pas_read_distribution.pdf"
    c. A scatterplot comparing our data's PAS usage distribution to bulk 3P-seq data named "bulk_data_comparison_ks.pdf"
    d. A collection of inverted CDF plots of various randomly-selected genes (from 9 quantiles based on length and 
    expression level) named "bulk_data_comparison_icdf.pdf"
    e. A histogram of post-filtration gene expression levels named "gene_hist.pdf"
    f. A histogram of post-filtration cell read counts named "cell_hist.pdf"
    g. A hexbinned TSNE distribution of our cells colored by UTR deviations named "tsne_hex_utr.png"
    h. A hexbinned TSNE distribution of our cells colored by bin density named "tsne_hex_density.png" 
    i. A hexbinned TSNE distribution of our cells, separated by cell ages, colored by UTR deviations from gene mean UTR
    lengths within each age named "age_tsne_hex_age_utr.png"
    j. A hexbinned TSNE distribution of our cells, separated by cell ages, colored by UTR deviations from gene mean UTR
    lengths within the full dataset named "age_tsne_hex_age_utr.png"
    k. A hexbinned subcluster TSNE distribution of our cells, separated by cell clusters, colored by UTR deviations from 
    gene mean UTR lengths within each cluster named "sub_tsne_hex_cluster_utr.png"
    l. A hexbinned subcluster TSNE distribution of our cells, separated by cell clusters, colored by UTR deviations from 
    gene mean UTR lengths within the full dataset named "sub_tsne_hex_utr.png"
    m. A heatmap of gene utr length deviations by cell ages named "gene_age_heatmap.pdf"
    n. A heatmap of gene utr length deviations by cell clusters "gene_cluster_heatmap.pdf"
    o. A heatmap of gene utr length deviations within selected neuronal cell clusters named "neuron_dev_heatmap.pdf"
    p. A collection of ICDF plots of selected genes that have interesting or important PAS usage patterns over cell
    development. These plots are organized into subfolders in the 'icdf' folder based on designation of their relevance:
        pa. "lengthening" contains genes that lengthen significantly at each cell age.
        pb. "shortening" contains genes that shorten significantly at each cell age.
        pc. "lengthening then shortening" contains genes that lengthen, then shorten significantly during development.
        pd. "shortening then lengthening" contains genes that shorten, then lengthen significantly during development.
        pe. "neuron_dev" are genes that show significant PAS usage pattern differences between the developmental 
        trajectories of excitatory and inhibitory neurons.
'''


def choose_transcripts(trans_to_gene_path):
    """Attaches each transcript to a gene based on gencode transcript annotations and chooses the most expressed
    transcript of each gene to use as our representative transcript of that gene."""

    # Reads in the transcript expression level data.
    transcript_counts = {}
    with open('transcript_counts_overlapfiles.txt', 'rt') as transcript_count_file:
        for line in transcript_count_file:
            transcript_data = line.rstrip().split()
            transcript_counts[transcript_data[1]] = int(transcript_data[0])

    # Associates each transcript to the gene it originated from
    gene_transcripts = {}
    gene_shorthands = {}
    with open('gencode.vM25.annotation.gtf', 'rt') as gene_annotation_file:
        reader = csv.reader(gene_annotation_file, delimiter='\t')
        # Skips the header
        for _ in range(0, 5):
            next(reader)
        for row in reader:
            gene_info = row[8]
            # Takes only rows that list a transcript and a gene id for a protein-coding gene.
            if 'ENSMUSG' not in gene_info or 'ENSMUST' not in gene_info or \
                    'transcript_type "protein_coding"' not in gene_info:
                continue
            gene_id = re.search('ENSMUSG.*?"', gene_info).group(0)[:-1]
            transcript_id = re.search('ENSMUST.*?"', gene_info).group(0)[:-1]
            gene_shorthand = re.search('gene_name ".*?"', gene_info).group(0)[11:-1]
            if gene_id not in gene_transcripts:
                gene_transcripts[gene_id] = [transcript_id]
                gene_shorthands[gene_id] = gene_shorthand
            elif transcript_id not in gene_transcripts[gene_id]:
                gene_transcripts[gene_id].append(transcript_id)
            else:
                continue

    trans_to_gene_out = open(trans_to_gene_path, 'wt')
    approved_transcripts_out = open(TRANSCRIPT_FILE_PATH, 'wt')

    # Chooses the most expressed transcript for each gene.
    for gene in gene_transcripts:
        most_expressed_transcript = ""
        most_expressed_count = 0
        for transcript in gene_transcripts[gene]:
            if transcript not in transcript_counts:
                continue
            transcript_count = transcript_counts[transcript]
            if transcript_count > most_expressed_count:  # and transcript_count > 1000:
                most_expressed_transcript = transcript
                most_expressed_count = transcript_count
        if most_expressed_transcript:
            shorthand = gene_shorthands[gene]
            trans_to_gene_out.write(most_expressed_transcript + '\t' + gene + '\t' + shorthand + '\n')
            approved_transcripts_out.write(most_expressed_transcript + '\n')

    trans_to_gene_out.close()
    approved_transcripts_out.close()


def generate_cell_data(cell_data_path):
    """Generates a dictionary of cell data keyed to each cell's name for internal use."""

    print("Generating cell data...")

    # Point this at your cell annotation file
    anno_path = "cell_annotations.csv"

    # Reformats given cell annotations for better python use.
    cell_data_dict = {}
    cell_id = -1
    with open(anno_path, 'rt') as cell_data_in:
        reader = csv.reader(cell_data_in, delimiter=',')
        for row in reader:
            # Skips the header row.
            if cell_id == -1:
                cell_id = 0
                continue
            cell = row[0]
            # To normalize gene expression levels.
            exon_count = row[1]
            age = row[8]
            tsne1 = row[13]
            tsne2 = row[14]
            subcluster = row[15]
            stsne1 = row[16]
            stsne2 = row[17]
            cluster = row[22]
            # Saves a lot of trouble in later analysis.
            if cluster == "Chondrocytes & osteoblasts":
                cluster = "Chondrocytes and osteoblasts"
            cluster = cluster.replace(' ', '_')
            trajectory = row[23]
            trajectory = trajectory.replace(' ', '_')
            if trajectory == "Neural_crest_1":
                trajectory = "Neural_crest_PNS_neurons_trajectory"
            elif trajectory == "Neural_crest_2":
                trajectory = "Neural_crest_PNS_glia_trajectory"
            elif trajectory == "Neural_crest_3":
                trajectory = "Neural_crest_melanocytes_trajectory"
            umap1 = row[24]
            umap2 = row[25]
            umap3 = row[26]
            rtrajectory = row[27]
            rumap1 = row[28]
            rumap2 = row[29]
            rumap3 = row[30]
            strajectory = trajectory + "." + row[31]
            strajectory = strajectory.replace(' ', '_')
            sumap1 = row[32]
            sumap2 = row[33]
            pseudotime = row[35]
            embryo_id = row[5]
            extraction_date = row[7]

            cell_data_dict[cell] = [age, exon_count, cluster, tsne1, tsne2, subcluster, stsne1, stsne2,
                                    trajectory, umap1, umap2, umap3, rtrajectory, rumap1, rumap2, rumap3,
                                    strajectory, sumap1, sumap2, cell_id, pseudotime, embryo_id, extraction_date]
            # Data indices
            # Age: 0
            # Exon Count: 1
            # Cluster: 2
            # TSNE Coordinates: 3, 4
            # Subcluster: 5
            # Sub TSNE Coordinates: 6, 7
            # Trajectory: 8
            # UMAP Coordinates: 9, 10, 11
            # Refined Trajectory: 12
            # Refined UMAP Coordinates: 13, 14, 15
            # Subtrajectory: 16
            # Sub UMAP Coordinates: 17, 18
            # Cell ID: 19
            # Pseudotime: 20
            # Embryo ID and Extraction Date: 21, 22

            cell_id += 1

    with open(cell_data_path, 'wb') as cell_data_out:
        pkl.dump(cell_data_dict, cell_data_out)

    print("Cell data generated.")


def initiate_data_files(names_info_path):
    """Initiates the data file that will hold our final data with indexed cell and gene names."""

    with open("cell_data_dict.pkl", 'rb') as cell_data_file:
        cell_data_dict = pkl.load(cell_data_file)

    # Reads the assigned index of each cell in our cell data dict and outputs them in that order to a list.
    cell_names = [""] * len(list(cell_data_dict.keys()))
    for cell_name in cell_data_dict:
        cell_idx = cell_data_dict[cell_name][19]
        cell_names[cell_idx] = cell_name

    # Point this at wherever this file is. It's simply a list of transcripts and the genes they correspond to.
    gene_name_path = 'trans_to_gene.txt'
    gene_dict = {}
    gene_names = []
    with open(gene_name_path, 'rt') as genes_in:
        reader = csv.reader(genes_in, delimiter='\t')
        for row in reader:
            gene_dict[row[0]] = row[1]
    for transcript in APPROVED_TRANSCRIPTS:
        gene_names.append(gene_dict[transcript])

    # This file holds an array of gene and transcript names indexed to their numeric ids, for the reverse referencing.
    # We also put our cell name information here.
    with open(names_info_path, 'wb') as names_out:
        pkl.dump([cell_names, APPROVED_TRANSCRIPTS, gene_names], names_out)

    print("Data files initiated.")


def generate_gene_data_dict(gene_data_path):
    """Initiates a dictionary that will hold our gene information for this dataset with gene numeric ids."""

    print("Generating gene data dict...")

    # This file is a bed file sorted by loci and then gene names.
    # This is generated with >sort -t$'\t' -k1,1 -k4,4 -k2,2n <input.bed > mm10_3utr_sorted.bed
    utr_exon_file_path = "3utr_overlap.bed"

    # This will allow us to reference gene id by gene name for mapping utr length coordinates.
    gene_data_dict = {}
    for idx, gene_name in enumerate(APPROVED_TRANSCRIPTS):
        gene_data_dict[gene_name] = [idx, []]
    # After this, gene_data_dict now has one entry, which is the gene index, as well as an empty list for exons.

    gene_data_temp = {}
    with open(utr_exon_file_path, 'rt') as exon_file:
        reader = csv.reader(exon_file, delimiter='\t')
        for row in reader:
            gene = row[3]
            if gene not in gene_data_dict:
                continue
            start = int(row[1])
            stop = int(row[2])
            length = stop - start
            if gene not in gene_data_temp:
                strand = row[5]
                gene_data_temp[gene] = [strand, []]
            gene_data_temp[gene][1].append([start, stop, length])
    # After this, our gene_data_temp is filled with strand and exon info of each gene.

    # Sorts our exons by ascending locus or descending locus for positive and negative stranded genes, respectively.
    # The "length" coordinate of each exon is now replaced with total length of preceding exons, for late calculations.
    for gene_temp in gene_data_temp:
        total_length = 0
        gene_exons = gene_data_temp[gene_temp][1]
        exons_sorted = []
        if gene_data_temp[gene_temp][0] == "+":
            for exon in sorted(gene_exons, key=lambda x: x[0]):
                exon_sorted = exon[0:2] + [total_length]
                total_length = total_length + exon[2]
                exons_sorted.append(exon_sorted)
        else:
            for exon in reversed(sorted(gene_exons, key=lambda x: x[0])):
                exon_sorted = exon[0:2] + [total_length]
                total_length = total_length + exon[2]
                exons_sorted.append(exon_sorted)
        gene_data_dict[gene_temp][1] = exons_sorted

    # This pickle file holds a dictionary which can be searched by gene name which will return the gene's numeric id.
    # This is done because searching the list of gene names for the target index is considerably slower.
    # We don't need such a file for our cells, as each cell already has this data attached to it in cell_age_dict.
    with open(gene_data_path, 'wb') as gene_data_out:
        pkl.dump(gene_data_dict, gene_data_out)

    print("Gene data dict generated.")


def parallel_job_submission(protocol):
    """Will submit an analysis of analysis_type to a folder named output_folder.
    The job script name is irrelevant, but set by this program internally for consistency."""

    sub_script = SCRIPT_PATH + protocol + "_matrixlengthandisoformanalysis.py"
    job_script = SCRIPT_PATH + "job_" + protocol + "_" + PAS_DATASET + ".sh"
    output_folder = PAS_DATASET + "/" + protocol + "/"
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    print("Generating job file...")
    # This bash script submits a parallel array of jobs that each run our subscript on one file in an input folder.
    # The subscript is called with two arguments: the path to its target file and the PAS dataset to use.
    with open(job_script, 'wt') as job_file:
        lines = ["#!/bin/bash\n",
                 "# This is a subfile of matrixlengthandisoformanalysis.py.\n",
                 "# This file is a job submission protocol for " + sub_script + ".\n",
                 "num=${SGE_TASK_ID}\n",
                 'printf -v num "%03d" $num\n',
                 '~sereno/software/anaconda2/bin/python ' + sub_script + ' ' + OVERLAP_PATH +
                 'utr_overlaps_${num}.bed.gz ' + PAS_DATASET + ' ' + os.getcwd() + '/\n']
        job_file.writelines(lines)
    print("Job file generated.")

    overlap_count = len(glob2.glob(OVERLAP_PATH + '*.bed.gz'))
    if len(os.listdir(output_folder)) == 0:
        print("Submitting job array...")
        # Limits the limit of jobs submitted at once to 10 if the total number of jobs is more than 10.
        # If the job is TSV writing or gene-by-cell calculations, raises limit to 20 to speed up process
        # (these jobs don't need much memory but will run fairly slowly due to high iteration counts).
        # These limits can be increased on clusters with more available memory (expect each job to use up to 20GB).
        if protocol == 'tsv' or protocol == 'raw' or protocol == 'centered':
            parallel_limit = str(30 if overlap_count >= 30 else overlap_count)
        else:
            parallel_limit = str(20 if overlap_count >= 20 else overlap_count)
        job_submission = "qsub -t 1-" + str(overlap_count) + " -tc " + parallel_limit + \
                         " -N " + protocol + " bash " + job_script
        print("Parallel submission protocol:")
        print(job_submission)
        process = subprocess.Popen(job_submission.split(), stdout=subprocess.PIPE)
        # These are required for correct subprocess functions in parallel.
        output, error = process.communicate()
        if output:
            print("Job submission report:")
            print(output.decode())
        if error:
            print("Job submission error:")
            print(error.decode())

        # Waits until we have as many output files as input files.
        file_check_flag = True
        while file_check_flag:
            if len(os.listdir(output_folder)) == overlap_count:
                file_check_flag = False
            time.sleep(300)
        print("All jobs finished.")
    elif len(os.listdir(output_folder)) != overlap_count:
        print("Something is wrong with your " + output_folder + " folder.")
        print("Remove all files from the folder and re-run this script.")
        sys.exit()
    else:
        print("Parallel processing output found.")


def generate_pas_files():
    """Builds a dictionary that holds info on each polyA site for each gene, as well as their relative read distance.
    Function now deprecated, code uses generate_pas_files_2 in its place. Left in for now."""

    print("Generating PAS data dict...")

    # If we provide a PAS dataset, this reads each entry in the data file and creates a dictionary of each annotated
    # PAS attached to each gene. This will be referenced later by our analysis to attach each read to an annotated PAS.
    pas_path = PAS_DATASET + "/"
    pas_file_path = pas_path + "PAS.bed"
    if not os.path.exists(pas_file_path):
        exit("PAS dataset missing. Ensure you have the file 'pas.bed' in the your " + pas_path + " folder.")

    with open("gene_data_dict.pkl", 'rb') as gene_data_in:
        gene_data_dict = pkl.load(gene_data_in)

    pas_data_dict = {}
    with open(pas_file_path, 'rt') as pas_file_in:
        reader = csv.reader(pas_file_in, delimiter='\t')
        for row in reader:
            gene = row[9]
            if gene not in APPROVED_TRANSCRIPTS:
                continue

            chrom = row[0]
            dataset = row[3]
            strand = row[5]
            if strand == "+":
                locus = int(row[2])
            else:
                locus = int(row[1])

            if chrom not in pas_data_dict:
                pas_data_dict[chrom] = {}
            if strand not in pas_data_dict[chrom]:
                pas_data_dict[chrom][strand] = {}
            if gene not in pas_data_dict[chrom][strand]:
                pas_data_dict[chrom][strand][gene] = []
            # Note that the length calculation doesn't need to take into account intronic regions, as the PAS data
            # file already reports which specific 3' UTR exon each PAS overlaps.
            if strand == "+":
                gene_locus = int(row[7])
                iso_length = locus - gene_locus
            else:
                gene_locus = int(row[8])
                iso_length = gene_locus - locus
            # Adds the total length of all previous exons to our PAS length.
            pre_exon_length = 0
            for exon in gene_data_dict[gene][1]:
                no_exon_overlap = True
                pre_exon_length = pre_exon_length + exon[2]
                if exon[0] <= locus <= exon[1]:
                    total_length = iso_length + pre_exon_length
                    no_exon_overlap = False
                    break
            # If no 3'UTR exon overlaps are found, the PAS is discarded.
            if no_exon_overlap:
                continue

            pas_data_dict[chrom][strand][gene].append([locus, total_length, dataset])

    # This is used as-is the function mapping script. It will be modified after function mapping is completed.
    with open(pas_path + "pas_data_dict.pkl", 'wb') as pas_file_out:
        pkl.dump(pas_data_dict, pas_file_out)

    # Adds our generated PAS to our name dictionary so we can later reference individual PAS by global indices.
    with open("names_by_id.pkl", 'rb') as names_in:
        pas_names = []
        names = pkl.load(names_in)
        for chrom in pas_data_dict:
            for strand in pas_data_dict[chrom]:
                for gene in pas_data_dict[chrom][strand]:
                    for pas_idx, pas in enumerate(pas_data_dict[chrom][strand][gene]):
                        pas_name = gene + "." + str(pas_idx)
                        pas_names.append(pas_name)
        names.append(pas_names)

    with open("names_by_id.pkl", 'wb') as names_out:
        pkl.dump(names, names_out)

    print("PAS data dict generated.")

    print("Generating PAS overlap function...")

    function_file_path = pas_path + "function/"
    parallel_job_submission("function")
    function_array = np.array([0] * 161)
    for function_slice_file_path in glob2.glob(function_file_path + "*.pkl"):
        with open(function_slice_file_path, 'rb') as function_file:
            function_slice = pkl.load(function_file)
            function_array += function_slice

    # Saves the full PAS read distribution we use for visualization
    pas_function = function_array / np.sum(function_array)
    function_output_path = pas_path + "pas_function_full.pkl"
    with open(function_output_path, 'wb') as function_output:
        pkl.dump(pas_function, function_output)

    function_axis = range(-400, 401, 5)
    plt.plot(function_axis, pas_function)
    plt.grid()
    plt.xlabel('Distance from PAS')
    plt.ylabel('Read Proportion')
    plt.title('Read Distribution in Relation to Annotated PAS')
    plt.savefig(pas_path + "figures/pas_read_distribution.pdf")
    plt.close()

    # Outputs the PAS read distribution we reference for assigning reads to PAS
    # This subset will map the indices [0, 161] to [+25, -201] for later calculations.
    function_implement = function_array[40:86]
    function_implement_path = pas_path + "pas_function.pkl"
    with open(function_implement_path, 'wb') as function_implement_out:
        pkl.dump(function_implement[::-1], function_implement_out)

    print("PAS read distribution plot generated.")


def generate_pas_files_2():
    """Builds a dictionary that holds info on each polyA site for each gene, as well as their relative read distance."""
    # Placeholder function for new distribution of reads.

    print("Generating PAS data dict...")

    # If we provide a PAS dataset, this reads each entry in the data file and creates a dictionary of each annotated
    # PAS attached to each gene. This will be referenced later by our analysis to attach each read to an annotated PAS.
    pas_path = PAS_DATASET + "/"
    pas_file_path = pas_path + "PAS.bed"
    if not os.path.exists(pas_file_path):
        exit("PAS dataset missing. Ensure you have the file 'pas.bed' in the your " + pas_path + " folder.")

    with open("gene_data_dict.pkl", 'rb') as gene_data_in:
        gene_data_dict = pkl.load(gene_data_in)

    pas_data_dict = {}
    with open(pas_file_path, 'rt') as pas_file_in:
        reader = csv.reader(pas_file_in, delimiter='\t')
        for row in reader:
            gene = row[9]
            if gene not in APPROVED_TRANSCRIPTS:
                continue

            chrom = row[0]
            dataset = row[3]
            strand = row[5]
            if strand == "+":
                locus = int(row[2])
            else:
                locus = int(row[1])

            if chrom not in pas_data_dict:
                pas_data_dict[chrom] = {}
            if strand not in pas_data_dict[chrom]:
                pas_data_dict[chrom][strand] = {}
            if gene not in pas_data_dict[chrom][strand]:
                pas_data_dict[chrom][strand][gene] = []
            # Note that the length calculation doesn't need to take into account intronic regions, as the PAS data
            # file already reports which specific 3' UTR exon each PAS overlaps.
            if strand == "+":
                gene_locus = int(row[7])
                iso_length = locus - gene_locus
            else:
                gene_locus = int(row[8])
                iso_length = gene_locus - locus
            # Adds the total length of all previous exons to our PAS length.
            pre_exon_length = 0
            for exon in gene_data_dict[gene][1]:
                no_exon_overlap = True
                pre_exon_length = pre_exon_length + exon[2]
                if exon[0] <= locus <= exon[1]:
                    total_length = iso_length + pre_exon_length
                    no_exon_overlap = False
                    break
            # If no 3'UTR exon overlaps are found, the PAS is discarded.
            if no_exon_overlap:
                continue

            pas_data_dict[chrom][strand][gene].append([locus, total_length, dataset, 0])

    # This is used as-is the function mapping script. It will be modified after function mapping is completed.
    with open(pas_path + "pas_data_dict.pkl", 'wb') as pas_file_out:
        pkl.dump(pas_data_dict, pas_file_out)

    # Adds our generated PAS to our name dictionary so we can later reference individual PAS by global indices.
    with open("names_by_id.pkl", 'rb') as names_in:
        pas_names = []
        names = pkl.load(names_in)
        for chrom in pas_data_dict:
            for strand in pas_data_dict[chrom]:
                for gene in pas_data_dict[chrom][strand]:
                    for pas_idx, pas in enumerate(pas_data_dict[chrom][strand][gene]):
                        pas_name = gene + "." + str(pas_idx)
                        pas_names.append(pas_name)
        names.append(pas_names)

    with open("names_by_id.pkl", 'wb') as names_out:
        pkl.dump(names, names_out)

    print("Gene name data generated.")

    print("Generating PAS overlap function...")

    # Count the number of reads within allowed range of each PAS.
    assignment_file_path = pas_path + "assignment/"
    parallel_job_submission("assignment")
    assignment_file_count = 0
    for assignment_slice_file_path in glob2.glob(assignment_file_path + "*.pkl"):
        with open(assignment_slice_file_path, 'rb') as assignment_file:
            pas_data_slice = pkl.load(assignment_file)
            if assignment_file_count == 0:
                pas_data_dict = pas_data_slice
            else:
                for chrom in pas_data_dict:
                    for strand in pas_data_dict[chrom]:
                        for gene in pas_data_dict[chrom][strand]:
                            for pas_idx, pas_data in enumerate(pas_data_dict[chrom][strand][gene]):
                                pas_count = pas_data_dict[chrom][strand][gene][pas_idx][3]
                                # Adds this slice's counts to the head dictionary
                                pas_count += pas_data_slice[chrom][strand][gene][pas_idx][3]
                                pas_data_dict[chrom][strand][gene][pas_idx][3] = pas_count
        assignment_file_count += 1
        print("Assignment file " + str(assignment_file_count) + " processed.")

    with open(pas_path + "pas_data_dict.pkl", 'wb') as pas_file_out:
        pkl.dump(pas_data_dict, pas_file_out)

    print("PAS data dict generated.")

    # Note for testing purposes that re-running of the assignment protocol is not necessary for distribution mapping.
    # However, it is necessary if you wish to build a mapping function from the distribution mapping, or the mapping
    # only allows one entry from each read.
    function_file_path = pas_path + "distribution/"
    parallel_job_submission("distribution")
    # Use 161 as the length for -400 to 400, 86 for -400 to 25 and 65 for -300 to 20
    # Remember to change the values in the subscripts 'assignment', 'distribution', 'raw', and 'centered' as well.
    function_array = np.array([0] * 65)
    for function_slice_file_path in glob2.glob(function_file_path + "*.pkl"):
        with open(function_slice_file_path, 'rb') as function_file:
            function_slice = pkl.load(function_file)
            function_array += function_slice

    # Saves the full PAS read distribution we use for visualization
    pas_function = function_array / np.sum(function_array)
    function_output_path = pas_path + "pas_function_full.pkl"
    with open(function_output_path, 'wb') as function_output:
        pkl.dump(pas_function, function_output)

    # Use range of (min, max+1, 5)
    function_axis = range(-300, 21, 5)
    plt.plot(function_axis, pas_function)
    plt.grid()
    plt.xlabel('Distance from PAS')
    plt.ylabel('Read Proportion')
    plt.title('Read Distribution in Relation to Annotated PAS')
    plt.savefig(pas_path + "figures/pas_read_distribution.pdf")
    plt.close()

    print("PAS read distribution plot generated.")


def generate_pas_venn(pas_venn_path):
    """Make a venn diagram of our PAS datasets' site overlap"""
    plt.figure(figsize=(10, 10))
    # These numbers are from overlapping our three PAS databases using bedtools intersect (-s -u)
    venn = venn3(subsets=(51284, 32026, 41632, 17217, 4095, 4239, 31041),
                 set_labels=('PolyADBv3', 'PolyASite2', 'GencodeM25'))
    for idx, subset in enumerate(venn.subset_labels):
        venn.subset_labels[idx].set_visible(False)
    plt.savefig(pas_venn_path)
    plt.close()


def analyze_gene_coverage():
    """Analyzes the relative coverage of each gene from our filtered and unfiltered data, compares it to the 3pseq bulk
    data, and plots the relationships as well as 18 random genes of varying lengths and expression levels."""

    # Merges all of our locus data and outputs it as a bedgraph.
    locus_file_count = 0
    locus_dict = {}
    for locus_file_path in glob2.glob(PAS_DATASET + "/pasloci/*.pkl"):
        with open(locus_file_path, 'rb') as locus_file:
            locus_dict_slice = pkl.load(locus_file)
            if locus_file_count == 0:
                locus_dict = locus_dict_slice
            else:
                for chrom in locus_dict_slice:
                    for locus in locus_dict_slice[chrom]:
                        if locus not in locus_dict[chrom]:
                            locus_dict[chrom][locus] = locus_dict_slice[chrom][locus]
                        else:
                            locus_dict[chrom][locus] += locus_dict_slice[chrom][locus]
        locus_file_count += 1
        print("Locus file " + str(locus_file_count) + " processed.")

    with open(PAS_DATASET + "/pas_locus_coverage.bg", 'wt') as locus_out_file:
        for chrom in locus_dict:
            for locus in locus_dict[chrom]:
                count = str(locus_dict[chrom][locus])
                stop = str(locus)
                start = str(locus - 1)
                locus_out_file.write('\t'.join([chrom, start, stop, count]) + '\n')

    with open('names_by_id.pkl', 'rb') as names_in:
        names = pkl.load(names_in)
        transcript_names = names[1]
        gene_names = names[2]

    # Merges all of the coverage dictionaries built in parallel by the 'raw' subscript
    coverage_path = PAS_DATASET + '/coverage/'
    coverage_dict = {}
    coverage_file_count = 0
    for coverage_file_path in glob2.glob(coverage_path + '*.pkl'):
        with open(coverage_file_path, 'rb') as coverage_file:
            coverage_subdict = pkl.load(coverage_file)
            if coverage_file_count == 0:
                coverage_dict = coverage_subdict
            else:
                for gene in coverage_subdict:
                    if gene not in coverage_dict:
                        coverage_dict[gene] = coverage_subdict[gene]
                    else:
                        for dataset in coverage_subdict[gene]:
                            if dataset not in coverage_dict[gene]:
                                coverage_dict[gene][dataset] = coverage_subdict[gene][dataset]
                            else:
                                for locus in coverage_subdict[gene][dataset]:
                                    if locus in coverage_dict[gene][dataset]:
                                        coverage_dict[gene][dataset][locus] += coverage_subdict[gene][dataset][locus]
                                    else:
                                        coverage_dict[gene][dataset][locus] = coverage_subdict[gene][dataset][locus]
        coverage_file_count += 1
        print("Coverage file " + str(coverage_file_count) + " processed.")

    print("Calculating K-S values...")
    # Calculates length and (filtered) expression level for each gene.
    expressions = []
    lengths = []
    unfiltered_mid_list = []
    filtered_mid_list = []
    gene_ordered_list = []
    mid_genes_list = []
    datasets = ['filtered', 'unfiltered', '3pseq']
    # inf_counts = 0
    for gene in coverage_dict:
        # Skips this gene if it doesn't have bulk data or filtered data. This should happen very rarely.
        if '3pseq' not in coverage_dict[gene] or 'filtered' not in coverage_dict[gene]:
            continue
        #
        # Calculates the longest recorded site in the gene and sets this as the gene range.
        longest = max([max(list(coverage_dict[gene][dataset].keys())) for dataset in datasets])
        coverage_dict[gene]['length'] = longest
        lengths.append(longest)
        #
        # Totals the counts of each PAS from our filtered data to get total gene expression
        expression = sum([coverage_dict[gene]['filtered'][locus] for locus in coverage_dict[gene]['filtered']])
        coverage_dict[gene]['expression'] = expression
        expressions.append(expression)
        #
        # So we can reference these genes later in the order we calculate their lengths and expressions.
        gene_ordered_list.append(gene)
        #
        raw_sites = list(coverage_dict[gene]['unfiltered'].keys())
        our_sites = list(coverage_dict[gene]['filtered'].keys())
        bulk_sites = list(coverage_dict[gene]['3pseq'].keys())
        # raw_dist = []
        # our_dist = []
        # bulk_dist = []
        # for site, coverage in coverage_dict[gene]['unfiltered'].items():
        #    raw_dist += [site] * coverage
        # for site, coverage in coverage_dict[gene]['filtered'].items():
        #    our_dist += [site] * coverage
        # for site, coverage in coverage_dict[gene]['3pseq'].items():
        #    bulk_dist += [site] * round(coverage)
        # Calculate a uniform binned distribution for the union of the datasets' bins.
        total_sites = sorted(list(set().union(raw_sites, our_sites, bulk_sites)))
        raw_dist = [0] * len(total_sites)
        our_dist = [0] * len(total_sites)
        bulk_dist = [0] * len(total_sites)
        for site_idx, site in enumerate(total_sites):
            if site in raw_sites:
                raw_dist[site_idx] = coverage_dict[gene]['unfiltered'][site]
            if site in our_sites:
                our_dist[site_idx] = coverage_dict[gene]['filtered'][site]
            if site in bulk_sites:
                bulk_dist[site_idx] = coverage_dict[gene]['3pseq'][site]
        #
        # Make the uniform binned distributions into cumulative density functions and perform tests.
        # "MID" stands for Mean Integral Difference.
        raw_cdf = np.cumsum(raw_dist) / sum(raw_dist)
        our_cdf = np.cumsum(our_dist) / sum(our_dist)
        bulk_cdf = np.cumsum(bulk_dist) / sum(bulk_dist)
        unfiltered_mid = np.sum(np.abs(raw_cdf - bulk_cdf)) / len(our_cdf)
        filtered_mid = np.sum(np.abs(our_cdf - bulk_cdf)) / len(our_cdf)
        # unfiltered_ks = -np.log10(sps.ks_2samp(raw_cdf, bulk_cdf[1])
        # filtered_ks = -np.log10(sps.ks_2samp(our_cdf, bulk_cdf)[1])
        # unfiltered_ks = -np.log10(logrank_test(raw_dist, bulk_dist).p_value)
        # filtered_ks = -np.log10(logrank_test(our_dist, bulk_dist).p_value)
        if unfiltered_mid == float('inf') or filtered_mid == float('inf'):
            # inf_counts += 1
            continue
        unfiltered_mid_list.append(unfiltered_mid)
        filtered_mid_list.append(filtered_mid)
        mid_genes_list.append(gene)
        coverage_dict[gene]['mid'] = [unfiltered_mid, filtered_mid]

    with open(PAS_DATASET + '/gene_coverage_dict.pkl', 'wb') as gene_coverage_dict_out:
        pkl.dump(coverage_dict, gene_coverage_dict_out)
    print("Coverage data processed.")

    plt.figure(figsize=(5, 5))
    plt.xlabel("Unfiltered Mean Integral Difference")
    plt.ylabel("Filtered Mean Integral Difference")
    plt.xlim(0, 0.5)
    plt.ylim(0, 0.5)
    ks_dataframe = pd.DataFrame(data={"Unfiltered Mean Integral Difference": unfiltered_mid_list,
                                      "Filtered Mean Integral Difference": filtered_mid_list})
    sns.kdeplot(data=ks_dataframe, shade=True)
    ks_figure_path = PAS_DATASET + "/figures/bulk_data_comparison_kde.pdf"
    plt.savefig(ks_figure_path, bbox_inches='tight')
    plt.close()

    # Splits our genes into 9 sections based on low, medium, and high expression and length.
    percentiles = [(100 / 3), (200 / 3)]
    expression_cutoffs = np.percentile(np.asarray(expressions), percentiles)
    length_cutoffs = np.percentile(np.asarray(lengths), percentiles)
    slice_array = np.empty((3, 3), dtype=object)
    for coords, _ in np.ndenumerate(slice_array):
        slice_array[coords] = []
    for idx, gene in enumerate(gene_ordered_list):
        expression_idx = sum([expressions[idx] > cutoff for cutoff in expression_cutoffs])
        length_idx = sum([lengths[idx] > cutoff for cutoff in length_cutoffs])
        slice_array[expression_idx, length_idx].append(gene)
    selected_genes = [[[] for _ in range(3)] for _ in range(3)]
    for coords, gene_list in np.ndenumerate(slice_array):
        selected_genes[coords[0]][coords[1]] = rd.sample(gene_list, 2)

    # Plots our data across 18 subplots.
    handles = []
    labels = []
    fig = plt.figure(figsize=(20, 20))
    outer = gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2)
    for i in range(9):
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i], wspace=0.1, hspace=0.1)
        for j in range(2):
            ax = plt.Subplot(fig, inner[j])
            gene = selected_genes[i // 3][i % 3][j]
            if i % 3 != 0:
                ax.axes.get_yaxis().set_visible(False)
            gene_data = coverage_dict[gene]
            length = gene_data['length']
            raw_sites = sorted(list(gene_data['unfiltered'].keys()))
            raw_entries = np.asarray([gene_data['unfiltered'][site] for site in raw_sites])
            our_sites = sorted(list(gene_data['filtered'].keys()))
            our_entries = np.asarray([gene_data['filtered'][site] for site in our_sites])
            bulk_sites = sorted(list(gene_data['3pseq'].keys()))
            bulk_entries = np.asarray([gene_data['3pseq'][site] for site in bulk_sites])
            # Makes an ICDF of our raw data.
            raw_sites = np.append(0, raw_sites)
            raw_sites = np.append(raw_sites, length + 1)
            raw_entries = raw_entries / sum(raw_entries)
            raw_cumulative = np.cumsum(raw_entries)
            raw_cumulative = np.append(np.asarray([0]), raw_cumulative)
            raw_function = [1 - n for n in raw_cumulative]
            raw_function = np.append(raw_function[0], raw_function)
            ax.step(raw_sites, raw_function, color='red', label='Our Data (Unfiltered)')
            # Makes an ICDF of our processed data
            our_sites = np.append(0, our_sites)
            our_sites = np.append(our_sites, length + 1)
            our_entries = our_entries / sum(our_entries)
            our_cumulative = np.cumsum(our_entries)
            our_cumulative = np.append(np.asarray([0]), our_cumulative)
            our_function = [1 - n for n in our_cumulative]
            our_function = np.append(our_function[0], our_function)
            ax.step(our_sites, our_function, color='green', label='Our Data (Filtered)')
            # Makes an ICDF of 3pseq data.
            bulk_sites = np.append(0, bulk_sites)
            bulk_sites = np.append(bulk_sites, length + 1)
            bulk_entries = bulk_entries / sum(bulk_entries)
            bulk_cumulative = np.cumsum(bulk_entries)
            bulk_cumulative = np.append(np.asarray([0]), bulk_cumulative)
            bulk_function = [1 - n for n in bulk_cumulative]
            bulk_function = np.append(bulk_function[0], bulk_function)
            ax.step(bulk_sites, bulk_function, color='blue', label='3Pseq Data')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            fig.add_subplot(ax)
            handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', frameon=False)
    icdf_figure_path = PAS_DATASET + "/figures/bulk_data_comparison_exp_age.pdf"
    plt.savefig(icdf_figure_path, bbox_inches='tight')
    plt.close()

    # Splits our genes into 9 sections based on low, medium, and high unfiltered and filtered MIDs.
    # We select these cutoffs to line up with the KDE density clouds.
    unfiltered_mid_cutoffs = [(.15, .2), (.34, .36),
                              (.1, .2), (.2, .3),
                              (.1, .15), (.15, .2), (.2, .25), (.25, .3), (.3, .35)]
    filtered_mid_cutoffs = [(.25, .3), (.25, .3),
                            (.15, .2), (.15, .2),
                            (0., .05), (0., .05), (0., .05), (0., .05), (0., .05)]
    cutoff_idxs = [i for i in range(9)]
    slice_array = [[] for _ in range(9)]
    # The mismatched boolean comparisons here are to better match the paneling to the scatter plot we generate.
    for gene, unfiltered_mid, filtered_mid in zip(mid_genes_list, unfiltered_mid_list, filtered_mid_list):
        for cutoff_idx, unfilt_co, filt_co in zip(cutoff_idxs, unfiltered_mid_cutoffs, filtered_mid_cutoffs):
            if unfilt_co[0] < unfiltered_mid < unfilt_co[1] and filt_co[0] < filtered_mid < filt_co[1]:
                slice_array[cutoff_idx].append(gene)
    selected_genes = ["" for _ in range(9)]
    for idx, gene_list in enumerate(slice_array):
        selected_genes[idx] = rd.sample(gene_list, 1)[0]

    selected_unfiltered_mid = []
    selected_filtered_mid = []
    handles = []
    labels = []
    fig, axs = plt.subplots(3, 3, sharey=True, figsize=(30, 15))
    for idx, gene in enumerate(selected_genes):
        row_idx = idx // 3
        col_idx = idx % 3
        gene_data = coverage_dict[gene]
        length = gene_data['length']
        gene_idx = mid_genes_list.index(gene)
        unfiltered_mid = unfiltered_mid_list[gene_idx]
        filtered_mid = filtered_mid_list[gene_idx]
        selected_unfiltered_mid.append(unfiltered_mid)
        selected_filtered_mid.append(filtered_mid)
        raw_sites = sorted(list(gene_data['unfiltered'].keys()))
        raw_entries = np.asarray([gene_data['unfiltered'][site] for site in raw_sites])
        our_sites = sorted(list(gene_data['filtered'].keys()))
        our_entries = np.asarray([gene_data['filtered'][site] for site in our_sites])
        bulk_sites = sorted(list(gene_data['3pseq'].keys()))
        bulk_entries = np.asarray([gene_data['3pseq'][site] for site in bulk_sites])
        # Makes an ICDF of our raw data.
        raw_sites = np.append(0, raw_sites)
        raw_sites = np.append(raw_sites, length + 1)
        raw_entries = raw_entries / sum(raw_entries)
        raw_cumulative = np.cumsum(raw_entries)
        raw_cumulative = np.append(np.asarray([0]), raw_cumulative)
        raw_function = [raw_cumulative[-1] - n for n in raw_cumulative]
        raw_function = np.append(raw_function[0], raw_function)
        axs[row_idx, col_idx].step(raw_sites, raw_function, color='red', label='Our Data (Unfiltered)')
        # Makes an ICDF of our processed data
        our_sites = np.append(0, our_sites)
        our_sites = np.append(our_sites, length + 1)
        our_entries = our_entries / sum(our_entries)
        our_cumulative = np.cumsum(our_entries)
        our_cumulative = np.append(np.asarray([0]), our_cumulative)
        our_function = [our_cumulative[-1] - n for n in our_cumulative]
        our_function = np.append(our_function[0], our_function)
        axs[row_idx, col_idx].step(our_sites, our_function, color='green', label='Our Data (Filtered)')
        # Makes an ICDF of 3pseq data.
        bulk_sites = np.append(0, bulk_sites)
        bulk_sites = np.append(bulk_sites, length + 1)
        bulk_entries = bulk_entries / sum(bulk_entries)
        bulk_cumulative = np.cumsum(bulk_entries)
        bulk_cumulative = np.append(np.asarray([0]), bulk_cumulative)
        bulk_function = [bulk_cumulative[-1] - n for n in bulk_cumulative]
        bulk_function = np.append(bulk_function[0], bulk_function)
        axs[row_idx, col_idx].step(bulk_sites, bulk_function, color='blue', label='3Pseq Data')
        gene_name = gene_names[transcript_names.index(gene)]
        axs[row_idx, col_idx].set_title(str(idx) + ': ' + gene_name)
        handles, labels = axs[row_idx, col_idx].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    icdf_figure_path = PAS_DATASET + "/figures/bulk_data_comparison_icdf_mid.pdf"
    plt.savefig(icdf_figure_path, bbox_inches='tight')
    plt.close()

    # Plots a scatter plot of all of our unfiltered and filtered MIDs, highlighting our specific genes.
    plt.figure(figsize=(5, 5))
    plt.xlabel("Unfiltered Mean Integral Difference")
    plt.ylabel("Filtered Mean Integral Difference")
    plt.xlim(0, 0.5)
    plt.ylim(0, 0.5)
    plt.scatter(unfiltered_mid_list, filtered_mid_list, marker='.', s=1)
    # Makes our selected genes' markers larger and red.
    plt.scatter(selected_unfiltered_mid, selected_filtered_mid, marker='.', s=10, c='r')
    # Circles our selected genes' markers
    # plt.scatter(selected_unfiltered_mid, selected_filtered_mid, s=10, facecolors='none', edgecolors='black')
    for i in range(9):
        plt.annotate(str(i), xy=(selected_unfiltered_mid[i], selected_filtered_mid[i]))
    ks_figure_path = PAS_DATASET + "/figures/bulk_data_comparison_scatter.pdf"
    plt.savefig(ks_figure_path, bbox_inches='tight')
    plt.close()


def calculate_averages(gene_average_path):
    """Creates a dictionary that holds the gene utr length averages of our given dataset."""

    print("Generating gene mean utr length data...")

    gene_average_dict = {}
    gene_dict = {}

    input_folder = PAS_DATASET + "/raw/"
    output_path = gene_average_path

    # These files are dictionaries returned by our raw data analysis subscript. Each one has two entries: a list of
    # cell, gene utr lengths and a subdictionary of each cell, gene utr length belonging to each subdivision of cells.
    # The subdivisions are cluster, trajectory, and age. Here we concatenate all of the data into one dictionary.
    file_count = 0
    for mean_file in glob2.glob(input_folder + "*.pkl"):
        with open(mean_file, "rb") as input_file:
            gene_subdict = pkl.load(input_file)
        for gene in gene_subdict:
            # If the gene hasn't been encountered yet, we simply append its data to our main dictionary.
            if gene not in gene_dict:
                gene_dict[gene] = gene_subdict[gene]
            # Otherwise, we concatenate once at a list level and once at a dictionary level
            else:
                gene_dict[gene][0] += gene_subdict[gene][0]
                for subdivision in gene_subdict[gene][1]:
                    if subdivision not in gene_dict[gene][1]:
                        gene_dict[gene][1][subdivision] = gene_subdict[gene][1][subdivision]
                    else:
                        gene_dict[gene][1][subdivision] += gene_subdict[gene][1][subdivision]
        file_count += 1
        print("Raw UTR Data File " + str(file_count) + " Processed.")

    # This calculates the average of each vector of information in the gene dict and outputs them to a mirror of the
    # dictionary, but with each key corresponding to averages across the original dataset rather than the dataset.
    for gene in gene_dict:
        gene_average_dict[gene] = [np.mean(gene_dict[gene][0]), {}]
        for subdivision in gene_dict[gene][1]:
            gene_average_dict[gene][1][subdivision] = np.mean(gene_dict[gene][1][subdivision])

    # Stores our gene average dictionary to be referenced by our data centering script.
    with open(output_path, 'wb') as gene_average_out:
        pkl.dump(gene_average_dict, gene_average_out)

    print("Gene mean utr length data generated.")


def analyze_pas_usage():
    """Analyzes the relative usage of each PAS by each cell, and formats this data for differential expression analysis.
    The differential expression analysis is performed using monocle in an external R script."""

    # Notes down highly expressed transcripts, for our purposes we use transcripts that have more than 10,000 counts in
    # our original overlap files (so 10,000 reads that overlap a 3'UTR).
    highly_expressed_transcripts_out = open('highly_expressed_transcripts.txt', 'wt')
    with open('transcript_counts_overlapfiles.txt', 'rt') as transcript_counts_in:
        for line in transcript_counts_in:
            transcript = line.split()[1]
            count = int(line.split()[0])
            if count >= 10000:
                highly_expressed_transcripts_out.write(transcript + '\n')
    highly_expressed_transcripts_out.close()

    parallel_job_submission("pasmatrix")
    
    exit("Don't forget to remove this exit statement.")

    # Loads in our cell data so we have a reference for which cells exist in the dataset and writes our cell metadata in
    # a way that will work for monocle3 differential analysis of certain cell features.
    cell_metadata_path = PAS_DATASET + "/cell_metadata.tsv"
    if os.path.exists(cell_metadata_path):
        print("Cell metadata found.")
    else:
        print("Generating cell metadata...")
        with open('cell_data_dict.pkl', 'rb') as cell_data_in:
            cell_data_dict = pkl.load(cell_data_in)
        with open('names_by_id.pkl', 'rb') as names_in:
            cell_names = pkl.load(names_in)[0]
        with open(cell_metadata_path, 'wt') as cell_metadata_out:
            cell_metadata_out.write("cell\tage\tcluster\ttrajectory\tbatch\tembryo\n")
            for cell in cell_names:
                cell_data = cell_data_dict[cell]
                cell_age = cell_data[0]
                cell_cluster = cell_data[2]
                cell_trajectory = cell_data[8]
                cell_embryo = cell_data[21]
                cell_extraction_date = cell_data[22]
                cell_metadata_out.write('\t'.join([cell, cell_age, cell_cluster,
                                                   cell_trajectory, cell_embryo, cell_extraction_date]) + '\n')
        print("Cell metadata generated.")

    with open('highly_expressed_transcripts.txt', 'rt') as highly_expressed_transcripts_in:
        highly_expressed_transcripts = [line.rstrip('\n') for line in highly_expressed_transcripts_in]
    # selected_transcripts = rd.sample(highly_expressed_transcripts, 100)
    selected_transcripts = highly_expressed_transcripts

    # Makes a hdf5 file that will store our matrices.

    # Reads through each pas matrix file and concatenates data for each gene. Creates pandas dataframes of each gene
    # in this way. Each data slice has unique cells so we don't need to check for existing cells each time.
    pas_matrix = {}
    pas_matrix_file_count = 0
    for pas_matrix_file_path in glob2.glob(PAS_DATASET + '/pasmatrix/*.pkl'):
        with open(pas_matrix_file_path, 'rb') as pas_matrix_file:
            pas_matrix_slice = pkl.load(pas_matrix_file)
            for transcript in selected_transcripts:
                if transcript not in pas_matrix_slice:
                    continue
                if transcript not in pas_matrix:
                    pas_matrix[transcript] = {}
                transcript_matrix = pas_matrix_slice[transcript]
                for coords in transcript_matrix:
                    pas_matrix[transcript][coords] = transcript_matrix[coords]
        pas_matrix_file_count += 1
        print("Pas matrix file " + str(pas_matrix_file_count) + " processed.")

    # Builds a dataframe using the lists of cell-pas expression pairs found in our main data. This is added to our
    # hdf5 file to be read by R, named after the transcript that sourced the data.
    gene_matrix_out_path = PAS_DATASET + '/pas_matrices.h5'
    gene_matrices_out = h5.File(gene_matrix_out_path, 'a')
    transcript_count = 0
    pas_info = []
    for transcript in pas_matrix:
        # If we've analyzed this transcript before, no need to do it again.
        if transcript in gene_matrices_out:
            continue
        transcript_matrix = pas_matrix[transcript]
        # Filters out transcripts that didn't have many PAS hits.
        if len(transcript_matrix) < 10000:
            continue
        cell_idxs = []
        pas_idxs = []
        expressions = []
        for coords in transcript_matrix:
            if coords == 'pas':
                pas_info = [pas_length for pas_length in transcript_matrix[coords]]
                continue
            cell_idx = coords[0]
            cell_idxs.append(cell_idx)
            pas_idx = coords[1]
            pas_idxs.append(pas_idx)
            expression = transcript_matrix[coords]
            expressions.append(expression)
        gene_matrices_out.create_dataset(transcript + "/cells", data=cell_idxs)
        gene_matrices_out.create_dataset(transcript + "/pas", data=pas_idxs)
        gene_matrices_out.create_dataset(transcript + "/exp", data=expressions)
        gene_matrices_out.create_dataset(transcript + "/pas.info", data=pas_info)
        transcript_count += 1
        if transcript_count % 100 == 0:
            print(str(transcript_count) + " transcripts processed.")
    gene_matrices_out.close()


def differential_pas_test():
    """"""

    # Re-orients our PAS data to fit with the new analysis. The other structure is faster for other analyses so we
    # maintain it
    refactored_pas_dict_path = PAS_DATASET + '/pas_data_by_tx.pkl'
    if os.path.exists(refactored_pas_dict_path):
        with open(refactored_pas_dict_path, 'rb') as pas_file_in:
            pas_data_by_transcript = pkl.load(pas_file_in)
    else:
        pas_data_by_transcript = {}
        with open(PAS_DATASET + '/pas_data_dict.pkl', 'rb') as pas_file_in:
            pas_data_dict = pkl.load(pas_file_in)
        for chrom in pas_data_dict:
            for strand in pas_data_dict[chrom]:
                for transcript in pas_data_dict[chrom][strand]:
                    pas_data_by_transcript[transcript] = pas_data_dict[chrom][strand][transcript]
        with open(refactored_pas_dict_path, 'wb') as pas_file_out:
            pkl.dump(pas_data_by_transcript, pas_file_out)

    transcript_frame_dict = {}
    for transcript in pas_data_by_transcript:
        transcript_frame_dict[transcript] = {}
        transcript_frame_dict[transcript]['ages'] = np.zeros((len(pas_data_by_transcript[transcript]), 5))
        transcript_frame_dict[transcript]['clusters'] = np.zeros((len(pas_data_by_transcript[transcript]), 38))

    with open("names_by_id.pkl", 'rb') as names_in:
        cell_names = pkl.load(names_in)[0]
    with open("cell_data_dict.pkl", 'rb') as cell_data_in:
        cell_data_dict = pkl.load(cell_data_in)

    cell_categories = {}
    for cell_idx, cell_name in enumerate(cell_names):
        cell_data = cell_data_dict[cell_name]
        cell_age = cell_data[0]
        if cell_age not in AGES:
            age_idx = 'NA'
        else:
            age_idx = AGES.index(cell_age)
        cell_cluster = cell_data[2]
        if cell_cluster not in CELL_TYPES:
            cluster_idx = 'NA'
        else:
            cluster_idx = CELL_TYPES.index(cell_cluster)
        cell_categories[cell_idx] = [age_idx, cluster_idx]

    pas_matrix_count = 0
    for pas_matrix_file_path in glob2.glob(PAS_DATASET + '/pasmatrix/*.pkl'):
        with open(pas_matrix_file_path, 'rb') as pas_matrix_file:
            pas_matrix = pkl.load(pas_matrix_file)
            for transcript in pas_matrix:
                transcript_data = pas_matrix[transcript]
                for coords in transcript_data:
                    # Skips data header.
                    if coords == 'pas':
                        continue
                    cell_idxs = cell_categories[coords[0]]
                    pas_idx = coords[1]
                    count = transcript_data[coords]
                    age_idx = cell_idxs[0]
                    cluster_idx = cell_idxs[1]
                    if age_idx == 'NA' or cluster_idx == 'NA':
                        continue
                    transcript_frame_dict[transcript]['ages'][(pas_idx, age_idx)] += count
                    transcript_frame_dict[transcript]['clusters'][(pas_idx, age_idx)] += count
        pas_matrix_count += 1
        print("PAS matrix " + str(pas_matrix_count) + " processed.")

    for transcript in transcript_frame_dict:
        age_frame = transcript_frame_dict[transcript]['ages']
        age_expr = np.sum(age_frame, axis=0) / np.sum(age_frame)
        pas_expr = np.sum(age_frame, axis=1)
        expected_age_frame = np.multiply.outer(pas_expr, age_expr)
        age_chs = sps.chisquare(age_frame.T, expected_age_frame.T)


def count_transcripts(transcript_count_path):
    """Sums our individual counts of each transcript that overlaps our PAS dataset for later calculations."""
    transcript_counts = {}
    file_count = 0
    for transcript_count_file_path in glob2.glob(PAS_DATASET + "/count/*.pkl"):
        with open(transcript_count_file_path, 'rb') as transcript_counts_in:
            transcript_counts_slice = pkl.load(transcript_counts_in)
            # Simple quick dictionary concatenation.
            if file_count == 0:
                transcript_counts = transcript_counts_slice
            else:
                for transcript in transcript_counts_slice:
                    if transcript in transcript_counts:
                        transcript_counts[transcript] += transcript_counts_slice[transcript]
                    else:
                        transcript_counts[transcript] = transcript_counts_slice[transcript]
            file_count += 1
            print("Transcript count file " + str(file_count) + " processed.")

    # Writes our transcript counts to our output text file.
    with open(transcript_count_path, 'wt') as transcript_counts_out:
        for transcript in transcript_counts:
            transcript_data_string = transcript + "\t" + str(transcript_counts[transcript]) + "\n"
            transcript_counts_out.write(transcript_data_string)


def format_expression_data(gene_expression_path):
    """Formats our gene expression data and count files for later plotting in R."""
    # Finds every possible transcript of every protein coding gene.
    gene_transcripts = {}
    with open('gencode.vM25.annotation.gtf', 'rt') as gene_annotation_file:
        reader = csv.reader(gene_annotation_file, delimiter='\t')
        # Skips the header
        for _ in range(0, 5):
            next(reader)
        for row in reader:
            gene_info = row[8]
            # Takes only rows that list a transcript and a gene id for a protein-coding gene.
            if 'ENSMUSG' not in gene_info or 'ENSMUST' not in gene_info or \
                    'transcript_type "protein_coding"' not in gene_info:
                continue
            gene_id = re.search('ENSMUSG.*?"', gene_info).group(0)[:-1]
            transcript_id = re.search('ENSMUST.*?"', gene_info).group(0)[:-1]
            if gene_id not in gene_transcripts:
                gene_transcripts[gene_id] = [transcript_id]
            elif transcript_id not in gene_transcripts[gene_id]:
                gene_transcripts[gene_id].append(transcript_id)
            else:
                continue

    # Reads in unfiltered transcript counts.
    transcript_counts_orig = {}
    with open('transcript_counts_original.txt', 'rt') as transcript_count_orig_file:
        for line in transcript_count_orig_file:
            transcript_data = line.rstrip().split()
            transcript_counts_orig[transcript_data[1]] = float(transcript_data[0])

    # Reads in filtered transcript counts.
    transcript_counts_pas = {}
    with open(PAS_DATASET + '/transcript_counts_pasoverlaps.txt', 'rt') as transcript_count_pas_file:
        for line in transcript_count_pas_file:
            transcript_data = line.rstrip().split()
            transcript_counts_pas[transcript_data[0]] = float(transcript_data[1])

    # Reads in baseline expression.
    gene_expr = {}
    with open('mouse_median_expr.txt', 'rt') as baseline_exp_file:
        for line in baseline_exp_file:
            gene_data = line.rstrip().split()
            gene_expr[gene_data[0]] = float(gene_data[1])

    gene_data = {}
    # This will iterate over transcripts in gencode25, our target reference set for coding genes.
    for gene in gene_transcripts:
        gene_base = gene.split('.')[0]
        # No baseline means we can't analyze this gene.
        if gene_base not in gene_expr:
            continue
        else:
            # Genes that we have no reads from but are in the gencode set default to 0 counts
            orig_expr = 0
            pas_expr = 0
            # Collates all of the transcripts of this gene isoform.
            for transcript in gene_transcripts[gene]:
                if transcript in transcript_counts_orig:
                    orig_expr += transcript_counts_orig[transcript]
                if transcript in transcript_counts_pas:
                    pas_expr += transcript_counts_pas[transcript]
            # If this is another isoform of an existing gene, we add its counts.
            # We don't add to the expression as this would double-count the baseline expression.
            if gene_base in gene_data:
                gene_data[gene_base] += np.asarray([0, orig_expr, pas_expr])
            else:
                base_expr = gene_expr[gene_base]
                gene_data[gene_base] = np.asarray([base_expr, orig_expr, pas_expr])

    with open(gene_expression_path, 'wt') as exp_info_out:
        # Header line for column names
        exp_info_out.write("gene\tbaseline\toriginal\tpas\n")
        # Formats all of our data as a string, appends gene name, and writes it to our file.
        for gene in gene_data:
            gene_data_string = '\t'.join([gene] + [str(data) for data in gene_data[gene]]) + '\n'
            exp_info_out.write(gene_data_string)


def plot_cluster_heatmap(cluster_gene_array, c_gene_ids_cutoff, filename):
    """Plots a heatmap of clusters-by-genes with whatever reduced array we pass it."""
    c_gene_array_cutoff = cluster_gene_array[c_gene_ids_cutoff, ]

    for c_coordinates, c_gene_age_entry in np.ndenumerate(c_gene_array_cutoff):
        if c_gene_age_entry < -50.0:
            c_gene_array_cutoff[c_coordinates] = -50.0
        elif c_gene_age_entry > 50.0:
            c_gene_array_cutoff[c_coordinates] = 50.0

    c_heatmap_path = PAS_DATASET + "/figures/" + filename
    c_heatmap_frame = pd.DataFrame(c_gene_array_cutoff,
                                   columns=[entry.replace('_', ' ') for entry in CELL_TYPES]).transpose()
    r_dist = sch.distance.pdist(c_heatmap_frame, metric='correlation')
    r_linkage = sch.linkage(r_dist, method='average')
    c_dist = sch.distance.pdist(c_heatmap_frame.transpose(), metric='correlation')
    c_linkage = sch.linkage(c_dist, method='average')
    c_sns_plot = sns.clustermap(c_heatmap_frame, cmap="viridis", xticklabels=False, figsize=(10, 10),
                                row_linkage=r_linkage, col_linkage=c_linkage)
    c_clustered_gene_ids = [c_gene_ids_cutoff[new_idx] for new_idx in c_sns_plot.dendrogram_col.reordered_ind]
    c_sns_plot.cax.set_visible(False)
    c_hm = c_sns_plot.ax_heatmap.get_position()
    c_sns_plot.ax_heatmap.set_position([c_hm.x0, c_hm.y0, c_hm.width * 4, c_hm.height])
    c_den = c_sns_plot.ax_col_dendrogram.get_position()
    c_sns_plot.ax_col_dendrogram.set_position([c_den.x0, c_den.y0, c_den.width * 4, c_den.height])
    c_sns_plot.savefig(c_heatmap_path)

    return c_clustered_gene_ids


def generate_icdfs(icdf_path, clustered_gene_ids, c_clustered_gene_ids_1, c_clustered_gene_ids_2):
    """Creates inverse cumulative density function plots of gene body coverage, representing relative usage of each PAS.
    Multiple ICDFs are plotted per figure, divided by some subdivision for individual genes with significant differences
    between PAS usages of each subdivision's cells."""

    print("Generating single-gene PAS usage plots...")

    minimum_exp = 1000
    if not os.path.exists(icdf_path):
        os.mkdir(icdf_path)

    with open("names_by_id.pkl", 'rb') as names_in:
        names = pkl.load(names_in)
        transcript_names = names[1]
    with open('trans_to_gene.txt', 'rt') as shorthand_in:
        reader = csv.reader(shorthand_in, delimiter='\t')
        short_names = []
        for row in reader:
            short_names.append(row[2])
        names.append(short_names)
    with open("names_by_id.pkl", 'wb') as names_out:
        pkl.dump(names, names_out)
    with open(PAS_DATASET + "/gene_age_array.pkl", 'rb') as arrays_in:
        arrays = pkl.load(arrays_in)
        age_transcript_array = arrays[0]
        age_count_array = arrays[1]
        cluster_gene_array = arrays[3]
        cluster_count_array = arrays[4]
    with open(PAS_DATASET + '/pas_data_dict.pkl', 'rb') as pas_data_in:
        pas_data_dict = pkl.load(pas_data_in)
    with open(PAS_DATASET + "/gene_average_dict.pkl", 'rb') as gene_averages_in:
        gene_average_dict = pkl.load(gene_averages_in)

    # Makes a list of the transcript ids whose transcripts experience an overall lengthening greater than
    # minimum_dev at each step of development, and have at least minimum_exp data points at each age.
    transcript_ids_total = []
    for transcript_id, transcript_row in enumerate(age_transcript_array):
        deviations = np.asarray([y - x for x, y in zip(transcript_row[:-1], transcript_row[1:])])
        if np.all(deviations >= 50) and np.all(age_count_array[transcript_id] >= minimum_exp):
            transcript_ids_total.append([transcript_id, "lengthening", AGES])
        if np.all(deviations <= -25) and np.all(age_count_array[transcript_id] >= minimum_exp):
            transcript_ids_total.append([transcript_id, "shortening", AGES])
        if (np.mean(deviations[2:4]) - np.mean(deviations[0:2])) >= 100 and \
                np.all(age_count_array[transcript_id] >= minimum_exp):
            transcript_ids_total.append([transcript_id, "shortening_then_lengthening", AGES])
        if (np.mean(deviations[2:4]) - np.mean(deviations[0:2])) <= -100 and \
                np.all(age_count_array[transcript_id] >= minimum_exp):
            transcript_ids_total.append([transcript_id, "lengthening_then_shortening", AGES])

    # Additionally, we plot all of the genes from our by-gene heatmaps.
    for clustered_gene_id in clustered_gene_ids:
        transcript_ids_total.append([clustered_gene_id, "age_dev_cov", AGES])

    for c_clustered_gene_id in c_clustered_gene_ids_1:
        transcript_ids_total.append([c_clustered_gene_id, "cluster_dev_1", CELL_TYPES])

    for c_clustered_gene_id in c_clustered_gene_ids_2:
        transcript_ids_total.append([c_clustered_gene_id, "cluster_dev_2", CELL_TYPES])

    # These transcript IDS correspond to Bdnf, Mmp9, Impa1, Braf1, Kcna1, Nrgn, Rhoa, Rspo3, Ascl1, Mapt, Calm1,
    # Ranbp1, Tcf4, and Nedd4l, respectively.
    for lit_id in [2171, 2822, 3042, 6745, 7442, 10943, 11664, 12002, 12509, 14205, 14946, 17741, 19715, 19683]:
        transcript_ids_total.append([lit_id, "lit_genes", CELL_TYPES])

    # These are indices of clusters that are involved in the developmental trajectory of specific neuron types.
    # Indices are Excitatory_neurons, Inhibitory_neuron_progenitors, Inhibitory_neurons, Neural_progenitor_cells,
    # Neural_Tube
    neuron_clusters = ['Cholinergic_neurons', 'Excitatory_neurons', 'Granule_neurons', 'Inhibitory_interneurons',
                       'Inhibitory_neuron_progenitors', 'Inhibitory_neurons', 'Neural_progenitor_cells', 'Neural_Tube',
                       'Postmitotic_premature_neurons', 'Sensory_neurons']
    neuron_dev_indices = [CELL_TYPES.index(cluster) for cluster in neuron_clusters]
    # neuron_dev = [CELL_TYPES[idx] for idx in [25, 24, 14, 10, 15]]
    neuron_gene_array = cluster_gene_array[:, neuron_dev_indices]
    neuron_array_means = np.mean(neuron_gene_array, axis=1)
    # Normalizes gene deviations to cluster expression levels.
    neuron_gene_array = (neuron_gene_array.T - neuron_array_means).T
    neuron_count_array = cluster_count_array[:, neuron_dev_indices]

    # Finds genes that contribute to differential lengthening between Excitatory and Inhibitory neurons.
    n_gene_ids_cutoff = []
    for n_gene_id, n_gene_row in enumerate(neuron_gene_array):
        if n_gene_id not in gene_average_dict:
            continue
        adjusted_mean = gene_average_dict[n_gene_id][0] + np.mean(neuron_gene_array)
        n_count_row = neuron_count_array[n_gene_id]
        standard_dev = np.std(n_gene_row)
        if standard_dev == 0:
            continue
        coefficient_of_variance = np.log10(standard_dev / adjusted_mean)
        if np.all(n_count_row >= 50) and coefficient_of_variance > -1.125:
            n_gene_ids_cutoff.append(n_gene_id)

    print(len(n_gene_ids_cutoff))

    n_gene_array_cutoff = neuron_gene_array[n_gene_ids_cutoff, ]

    # Old bugfix.
    # count = 0
    # n_gene_ids_cutoff_new = []
    # for id in n_gene_ids_cutoff:
    #     transcript = t_names_old[id]
    #     if transcript not in t_names_new:
    #         count += 1
    #     else:
    #         n_gene_ids_cutoff_new.append(id)

    # Constrains data to -100 to 100 range to eliminate outliers which will throw off analysis, then plots a heatmap.
    # This isn't done in the usual heatmap generation script as the arrays it references we use for icdf generation.
    for n_coordinates, n_gene_age_entry in np.ndenumerate(n_gene_array_cutoff):
        if n_gene_age_entry < -100.0:
            n_gene_array_cutoff[n_coordinates] = -100.0
        if n_gene_age_entry > 100.0:
            n_gene_array_cutoff[n_coordinates] = 100.0
    n_heatmap_path = PAS_DATASET + "/figures/neuron_dev_heatmap.pdf"
    n_heatmap_frame = pd.DataFrame(n_gene_array_cutoff,
                                   columns=[cluster.replace("_"," ") for cluster in neuron_clusters]).transpose()
    r_dist = sch.distance.pdist(n_heatmap_frame, metric='correlation')
    r_linkage = sch.linkage(r_dist, method='average')
    c_dist = sch.distance.pdist(n_heatmap_frame.transpose(), metric='correlation')
    c_linkage = sch.linkage(c_dist, method='average')
    n_sns_plot = sns.clustermap(n_heatmap_frame, cmap="viridis", xticklabels=False, figsize=(10, 10),
                                row_linkage=r_linkage, col_linkage=c_linkage)
    n_clustered_gene_ids = [n_gene_ids_cutoff[new_idx] for new_idx in n_sns_plot.dendrogram_col.reordered_ind]
    n_hm = n_sns_plot.ax_heatmap.get_position()
    n_sns_plot.ax_heatmap.set_position([n_hm.x0, n_hm.y0, n_hm.width * 2, n_hm.height])
    n_den = n_sns_plot.ax_col_dendrogram.get_position()
    n_sns_plot.ax_col_dendrogram.set_position([n_den.x0, n_den.y0, n_den.width * 2, n_den.height])
    n_sns_plot.savefig(n_heatmap_path)

    # Adds the gene ids important to neuron development to the transcript ids list.
    transcript_ids_total += [[transcript_id, "neuron_dev", neuron_clusters] for transcript_id in n_clustered_gene_ids]

    # Scans our PAS annotations for each transcript, then builds it a list of PAS site coordinates.
    transcript_site_dict = {}
    transcript_frame_dict = {}
    for transcript_id, transcript_data in enumerate(transcript_ids_total):
        transcript = transcript_names[transcript_data[0]]
        for chromosome in pas_data_dict:
            for strand in pas_data_dict[chromosome]:
                if transcript in pas_data_dict[chromosome][strand]:
                    if strand == "+":
                        transcript_site_dict[transcript_id] = \
                            [entry[1] for entry in pas_data_dict[chromosome][strand][transcript]]
                    else:
                        transcript_site_dict[transcript_id] = \
                            [entry[1] for entry in pas_data_dict[chromosome][strand][transcript]][::-1]
        # This builds an empty numpy array to put our data into, with rows corresponding to PAS isoforms
        # and columns corresponding to whatever subdivisions we're splitting the ICDF by.
        # We reference the transcript frame dict and site dict by internal indexing, but other objects by global
        # transcript/gene indexing. This is so we can build multiple frames across different axes for the same
        # transcript.
        transcript_frame_dict[transcript_id] = \
            np.zeros((len(transcript_site_dict[transcript_id]), len(transcript_data[2])))

    # Counts the occurrences of each PAS isoform at each subdivision, making an array for each transcript.
    file_count = 0
    for file_path in glob2.glob(PAS_DATASET + '/raw/*.pkl'):
        with open(file_path, 'rb') as pf:
            pas_usage_dict = pkl.load(pf)
            for transcript_id, transcript_data in enumerate(transcript_ids_total):
                if transcript_data[0] not in pas_usage_dict:
                    continue
                transcript_subdict = pas_usage_dict[transcript_data[0]][1]
                for sub_idx, sub in enumerate(transcript_data[2]):
                    if sub not in transcript_subdict:
                        continue
                    else:
                        for site in transcript_subdict[sub]:
                            transcript_sites = transcript_site_dict[transcript_id]
                            # If the PAS isoform isn't in our annotated sites (which means that this cell expressed
                            # the gene multiple times with different PAS usage, so we averaged the PAS coordinates),
                            # it attaches the data point to the nearest annotated site.
                            if site not in transcript_sites:
                                site = min(transcript_sites, key=lambda x: abs(x - site))
                            site_idx = transcript_sites.index(int(site))
                            transcript_frame_dict[transcript_id][(site_idx, sub_idx)] += 1
        file_count += 1
        print("PAS-by-age data file " + str(file_count) + " processed.")

    # viridis = cm.get_cmap('viridis')
    cluster_cm = plt.cm.get_cmap('nipy_spectral')
    # if not os.path.exists(icdf_path):
    #     os.mkdir(icdf_path)
    colors_age = ["#440154", "#3B528B", "#21918C", "#00A15B", "#00EC0A"]
    colors_neurons = [cluster_cm(x) for x in np.linspace(0, 1, 10)]
    colors_cluster = [cluster_cm(x) for x in np.linspace(0, 1, 38)]

    # Transcript_ids_total requires sublists of [transcript id, dataset designation, line designations].
    for transcript_id, transcript_data in enumerate(transcript_ids_total):
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        # Adds two extra sites: a '0' site to represent the start of our ICDF and an 'end' site to represent the end.
        sites = np.append(np.asarray([0]), np.asarray(transcript_site_dict[transcript_id]))
        sites = np.append(sites, np.asarray(sites[-1] + 1))
        transcript_designation = transcript_data[1]
        transcript_designation_directory = str(icdf_path) + transcript_designation + "/"
        # Makes a color map based on the number of subdivisions the icdf needs.
        # colors =
        if len(transcript_data[2]) == 5:
            colors = colors_age
        elif len(transcript_data[2]) == 10:
            colors = colors_neurons
        else:
            colors = colors_cluster
        if not os.path.exists(transcript_designation_directory):
            os.mkdir(transcript_designation_directory)
        for sub_idx, sub in enumerate(transcript_data[2]):
            entries = transcript_frame_dict[transcript_id][:, sub_idx]
            # Calculates the proportion of each PAS site, then makes a cumulative density function of them.
            entries = entries / sum(entries)
            cumulative = np.cumsum(entries)
            cumulative = np.append(np.asarray([0]), cumulative)
            # Inverts the cumulative density function
            function = [cumulative[-1] - n for n in cumulative]
            # Makes a flat starting segment between the '0' site and our first PAS.
            function = np.append(function[0], function)
            ax.step(sites, function, color=colors[sub_idx], label=sub.replace("_", " "))
        gene_name = short_names[transcript_data[0]]
        plt.title(gene_name + " PAS Usage")
        # Only plots legend if there is a manageable number of clusters.
        if len(transcript_data[2]) < 10:
            ax.legend(frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        if transcript_designation == "age_dev_cov":
            figure_path = transcript_designation_directory + str(clustered_gene_ids.index(transcript_data[0])) + "." + \
                          gene_name + ".png"
        elif transcript_designation == "neuron_dev":
            figure_path = transcript_designation_directory + str(n_clustered_gene_ids.index(transcript_data[0])) + "." \
                          + gene_name + ".png"
        elif transcript_designation == "cluster_dev_1":
            figure_path = transcript_designation_directory + str(c_clustered_gene_ids_1.index(transcript_data[0])) + \
                          "." + gene_name + ".png"
        elif transcript_designation == "cluster_dev_2":
            figure_path = transcript_designation_directory + str(c_clustered_gene_ids_2.index(transcript_data[0])) + \
                          "." + gene_name + ".png"
        else:
            figure_path = transcript_designation_directory + gene_name + ".png"
        fig.savefig(figure_path, bbox_inches='tight')
        plt.close('all')


def generate_heatmaps(heatmap_path):
    """Takes the numpy array pickle files stored in the heatmap/ subfolder of our PAS dataset folder created by our
    'heatmap' subscript and makes a single large pandas dataframe, using it to draw a gene by age heatmap. Also
    takes the count data stored in our h5 files and plots histograms of the data distribution."""

    print("Generating cell and gene expression statistics...")

    file_count = 0
    total_genes = len(APPROVED_TRANSCRIPTS)
    gene_counts = np.asarray([0] * total_genes)
    cell_counts = [0] * 2058652
    for centered_data_file in glob2.glob(PAS_DATASET + "/centered/*.h5"):
        with h5.File(centered_data_file, 'r') as count_input:
            gene_count_slice = np.asarray(list(count_input['gene_counts']))
            gene_counts += gene_count_slice
            # This is to avoid burdening our memory too heavily by repeatedly loading giant arrays.
            cell_count_data = list(count_input['cell_counts'])
            for cell_data in cell_count_data:
                cell_counts[cell_data[1]] = cell_data[0]
        file_count += 1
        print("Count File " + str(file_count) + " Processed.")

    cell_gene_expression_out = PAS_DATASET + "/cell_gene_expressions.pkl"
    with open(cell_gene_expression_out, 'wb') as cell_gene_expression_file:
        pkl.dump([cell_counts, gene_counts], cell_gene_expression_file)
    with open('names_by_id.pkl', 'rb') as names_in:
        names = pkl.load(names_in)

    pl.hist(gene_counts, bins=np.logspace(np.log10(1), np.log10(1000000), 50))
    pl.gca().set_xscale("log")
    pl.xlabel('Read Count')
    pl.ylabel('Gene Count')
    pl.title('Distributions of Read Counts per Gene')
    pl.savefig(PAS_DATASET + '/figures/gene_hist.pdf')
    pl.clf()
    pl.close()

    pl.hist(cell_counts, bins=np.logspace(np.log10(1), np.log10(2000), 50))
    pl.gca().set_xscale("log")
    pl.xlabel('Read Count')
    pl.ylabel('Cell Count')
    pl.title('Distributions of Read Counts per Cell')
    pl.savefig(PAS_DATASET + '/figures/cell_hist.pdf')
    pl.clf()
    pl.close()

    print("Cell and gene expression distributions plotted.")

    print("Generating gene by age heatmap...")

    file_count = 0
    for heatmap_file in glob2.glob(PAS_DATASET + "/heatmap/*.pkl"):
        with open(heatmap_file, 'rb') as array_input:
            arrays = pkl.load(array_input)
            gene_subarray = arrays[0]
            count_subarray = arrays[1]
            trajectory_subarray_dict = arrays[2]
            cluster_gene_subarray = arrays[3]
            cluster_count_subarray = arrays[4]
            # Initializes our gene by age array.
            if file_count == 0:
                gene_array = gene_subarray
                count_array = count_subarray
                trajectory_array_dict = trajectory_subarray_dict
                cluster_gene_array = cluster_gene_subarray
                cluster_count_array = cluster_count_subarray
            else:
                # The arrays are all generated to be the same length, and the coordinates will always match.
                for coordinates, old_mean in np.ndenumerate(gene_array):
                    old_count = count_array[coordinates]
                    new_mean = gene_subarray[coordinates]
                    new_count = count_subarray[coordinates]
                    total_count = old_count + new_count
                    if total_count == 0:
                        total_mean = 0
                    else:
                        total_mean = (old_mean * old_count + new_mean * new_count) / total_count
                    count_array[coordinates] = total_count
                    gene_array[coordinates] = total_mean
                for c_coordinates, c_old_mean in np.ndenumerate(cluster_gene_array):
                    c_old_count = cluster_count_array[c_coordinates]
                    c_new_mean = cluster_gene_subarray[c_coordinates]
                    c_new_count = cluster_count_subarray[c_coordinates]
                    c_total_count = c_old_count + c_new_count
                    if c_total_count == 0:
                        c_total_mean = 0
                    else:
                        c_total_mean = (c_old_mean * c_old_count + c_new_mean * c_new_count) / c_total_count
                    cluster_count_array[c_coordinates] = c_total_count
                    cluster_gene_array[c_coordinates] = c_total_mean
                for trajectory in trajectory_subarray_dict:
                    if trajectory not in trajectory_array_dict:
                        trajectory_array_dict[trajectory] = trajectory_subarray_dict[trajectory]
                        # This will probably never happen, but just in case a trajectory is entirely missing from our
                        # first cell set file.
                    else:
                        trajectory_gene_array = trajectory_array_dict[trajectory][0]
                        trajectory_count_array = trajectory_array_dict[trajectory][1]
                        trajectory_gene_subarray = trajectory_subarray_dict[trajectory][0]
                        trajectory_count_subarray = trajectory_subarray_dict[trajectory][1]
                        for coordinates, old_mean in np.ndenumerate(trajectory_gene_array):
                            old_count = trajectory_count_array[coordinates]
                            new_mean = trajectory_gene_subarray[coordinates]
                            new_count = trajectory_count_subarray[coordinates]
                            total_count = old_count + new_count
                            if total_count == 0:
                                total_mean = 0
                            else:
                                total_mean = (old_mean * old_count + new_mean * new_count) / total_count
                            trajectory_array_dict[trajectory][0][coordinates] = total_mean
                            trajectory_array_dict[trajectory][1][coordinates] = total_count
        file_count += 1
        print("Heatmap File " + str(file_count) + " Processed.")

    overframe_path = PAS_DATASET + "/gene_age_array.pkl"
    with open(overframe_path, 'wb') as overframe_out:
        pkl.dump([gene_array, count_array, trajectory_array_dict,
                  cluster_gene_array, cluster_count_array], overframe_out)

    with open(PAS_DATASET + "/gene_average_dict.pkl", 'rb') as gene_averages_in:
        gene_average_dict = pkl.load(gene_averages_in)

    gene_array_means = np.mean(gene_array, axis=1)
    # Normalizes gene deviations to age expression levels.
    gene_array = (gene_array.T - gene_array_means).T
    # Adjust this cutoff as it fits the data. We use 1000 as a basis, taking the upper 2/3 of our (pre-filtered) data.
    # We also take the top 1000 genes by standard deviation from these highly expressed genes, biased towards highly
    # expressed genes as these naturally will have lower variance.
    gene_ids_cutoff = []
    for gene_id, gene_row in enumerate(gene_array):
        if gene_id not in gene_average_dict:
            continue
        # Avoids divby0 errors. We wouldn't want such invariant genes in the heatmaps anyways.
        if np.std(gene_row) == 0:
            continue
        adjusted_mean = gene_average_dict[gene_id][0] + gene_array_means[gene_id]
        standard_dev = np.std(gene_array[gene_id])
        coefficient_of_variance = np.log10(standard_dev / adjusted_mean)
        # This cutoff selects the top 1000 genes by CoV
        if np.all(count_array[gene_id] >= 1000) and coefficient_of_variance > -1.4978:
            gene_ids_cutoff.append(gene_id)

    gene_array_cutoff = gene_array[gene_ids_cutoff, ]

    # Adjust these cutoffs as fits the data.
    for coordinates, gene_age_entry in np.ndenumerate(gene_array_cutoff):
        if gene_age_entry < -100.0:
            gene_array_cutoff[coordinates] = -100.0
        if gene_age_entry > 100.0:
            gene_array_cutoff[coordinates] = 100.0

    heatmap_frame = pd.DataFrame(gene_array_cutoff, columns=["9.5", "10.5", "11.5", "12.5", "13.5"])
    # Extracts the genes in the order they're clustered for our referencing.
    dist = sch.distance.pdist(heatmap_frame, metric='correlation')
    linkage = sch.linkage(dist, method='average')
    # Plots based off of our pre-computed clustering technique.
    sns_plot = sns.clustermap(heatmap_frame, yticklabels=False, col_cluster=False, cmap="viridis", figsize=(20, 20),
                              row_linkage=linkage)
    clustered_gene_ids = [gene_ids_cutoff[new_idx] for new_idx in sns_plot.dendrogram_row.reordered_ind]
    # Used to remove the legend for our manual figure editing.
    sns_plot.cax.set_visible(False)
    hm = sns_plot.ax_heatmap.get_position()
    sns_plot.ax_heatmap.set_position([hm.x0, hm.y0, hm.width * 0.25, hm.height])
    sns_plot.savefig(heatmap_path)

    cluster_array_means = np.mean(cluster_gene_array, axis=1)
    # Normalizes gene deviations to cluster expression levels.
    cluster_gene_array = (cluster_gene_array.T - cluster_array_means).T

    # Adjust all filters and cutoffs as they filter data. We used these cutoffs to filter for highly-expressed genes
    # that go through noticeable lengthening across different cell types.
    c_gene_ids_cutoff_1 = []
    c_gene_ids_cutoff_2 = []
    for c_gene_id, c_gene_row in enumerate(cluster_gene_array):
        if c_gene_id not in gene_average_dict:
            continue
        adjusted_mean = gene_average_dict[c_gene_id][0] + cluster_array_means[c_gene_id]
        c_count_row = cluster_count_array[c_gene_id]
        standard_dev = np.std(c_gene_row)
        if standard_dev == 0:
            continue
        coefficient_of_variance = np.log10(standard_dev / adjusted_mean)
        if np.all(c_count_row >= 50) and coefficient_of_variance > -1.696:
            c_gene_ids_cutoff_1.append(c_gene_id)
        if np.all(c_count_row >= 20) and coefficient_of_variance > -1.592:
            c_gene_ids_cutoff_2.append(c_gene_id)

    c_clustered_gene_ids_1 = plot_cluster_heatmap(cluster_gene_array, c_gene_ids_cutoff_1,
                                                  "cluster_gene_array_cutoff_50.pdf")
    c_clustered_gene_ids_2 = plot_cluster_heatmap(cluster_gene_array, c_gene_ids_cutoff_2,
                                                  "cluster_gene_array_cutoff_20.pdf")

    print("Heatmaps generated.")

    # Generates ICDFs from our interesting genes.
    icdf_path = PAS_DATASET + "/figures/icdfs/"

    if os.path.exists(icdf_path):
        print("ICDFs already generated.")
    else:
        print("Generating ICDFs...")
        generate_icdfs(icdf_path, clustered_gene_ids, c_clustered_gene_ids_1, c_clustered_gene_ids_2)


def calculate_transcript_stats(gene_stats_path):
    """Calculates the expression level of each transcript at each level of filtration: firstly, the raw bam data;
    secondly, the data that overlaps 3'UTR regions; and thirdly, the data that is assigned to annotated PAS."""

    print("Generating gene stats...")

    with open("names_by_id.pkl", 'rb') as name_file:
        names = pkl.load(name_file)
    transcripts = names[1]
    genes = names[2]

    total_genes = len(APPROVED_TRANSCRIPTS)
    gene_count_array = np.zeros((total_genes, 3))

    # First column: how many reads overlap any part of the entire transcript
    with open("transcript_counts_original.txt", 'rt') as transcript_file_original:
        for line in transcript_file_original:
            transcript = line.split()[1]
            if transcript not in transcripts:
                continue
            else:
                gene_id = transcripts.index(transcript)
                gene_count_array[gene_id, 0] = line.split()[0]

    # Second column: how many reads overlap the 3'UTR of the transcript
    with open("transcript_counts_overlapfiles.txt", 'rt') as transcript_file_original:
        for line in transcript_file_original:
            transcript = line.split()[1]
            if transcript not in transcripts:
                continue
            else:
                gene_id = transcripts.index(transcript)
                gene_count_array[gene_id, 1] = line.split()[0]

    # Third column: how many reads can be assigned an annotated PAS
    gene_age_array_path = PAS_DATASET + "/gene_age_array.pkl"
    with open(gene_age_array_path, 'rb') as gene_age_array_file:
        gene_age_array = pkl.load(gene_age_array_file)
        filtered_count_array = gene_age_array[1]
        filtered_count_vector = np.sum(filtered_count_array, axis=1)
        gene_count_array[:, 2] = filtered_count_vector

    gene_data_frame = pd.DataFrame(gene_count_array, index=genes)

    with open(gene_stats_path, 'w') as gene_data_out:
        gene_data_frame.to_csv(gene_data_out, sep='\t')

    print("Gene stats generated.")


def calculate_pas_stats(pas_stats_path):
    """Takes the sliced pas usage arrays generated by our 'centered' subscript and concatenates them, then outputs the
    concatenates array as tab-separated text files, labeled by PAS ids and clusters/ages."""

    print("Concatenating PAS usage stats...")

    with open("names_by_id.pkl", 'rb') as names_in:
        pas_names = pkl.load(names_in)[3]
        total_pas = len(pas_names)

    file_count = 0
    age_pas_array = np.zeros((total_pas, 5))
    cluster_pas_array = np.zeros((total_pas, 38))
    for pas_array_file_path in glob2.glob(PAS_DATASET + "/pasarrays/*.pkl"):
        with open(pas_array_file_path, 'rb') as pas_array_file:
            arrays = pkl.load(pas_array_file)
        if file_count == 0:
            age_pas_array = arrays[0]
            cluster_pas_array = arrays[1]
        else:
            age_pas_array += arrays[0]
            cluster_pas_array += arrays[1]
        file_count += 1
        print("PAS usage array file " + str(file_count) + " processed.")

    age_pas_frame = pd.DataFrame(age_pas_array, index=pas_names, columns=['9.5', '10.5', '11.5', '12.5', '13.5'])
    age_pas_out_path = PAS_DATASET + "/pas_stats_age.txt"
    with open(age_pas_out_path, 'w') as age_pas_out:
        age_pas_frame.to_csv(age_pas_out, sep='\t')

    cluster_pas_frame = pd.DataFrame(cluster_pas_array, index=pas_names,
                                     columns=['Cardiac_muscle_lineages', 'Cholinergic_neurons',
                                              'Chondroctye_progenitors', 'Chondrocytes_and_osteoblasts',
                                              'Connective_tissue_progenitors', 'Definitive_erythroid_lineage',
                                              'Early_mesenchyme', 'Endothelial_cells', 'Ependymal_cell',
                                              'Epithelial_cells', 'Excitatory_neurons', 'Granule_neurons',
                                              'Hepatocytes', 'Inhibitory_interneurons', 'Inhibitory_neuron_progenitors',
                                              'Inhibitory_neurons', 'Intermediate_Mesoderm', 'Isthmic_organizer_cells',
                                              'Jaw_and_tooth_progenitors', 'Lens', 'Limb_mesenchyme', 'Megakaryocytes',
                                              'Melanocytes', 'Myocytes', 'Neural_progenitor_cells', 'Neural_Tube',
                                              'Neutrophils', 'Notochord_cells', 'Oligodendrocyte_Progenitors',
                                              'Osteoblasts', 'Postmitotic_premature_neurons',
                                              'Premature_oligodendrocyte', 'Primitive_erythroid_lineage', 'Radial_glia',
                                              'Schwann_cell_precursor', 'Sensory_neurons', 'Stromal_cells',
                                              'White_blood_cells'])
    with open(pas_stats_path, 'w') as cluster_pas_out:
        cluster_pas_frame.to_csv(cluster_pas_out, sep='\t')


def generate_tsne_plots(tsne1, tsne2, utrs, ages, age_utrs, clusters, stsne1, stsne2, cluster_utrs):
    """Plots TSNE data based on our formatted single-cell data."""
    # This is the cutoff for our color bars, we use 50, near the 25th/75th percentile of our utr deviations, for this
    # analysis. Due to the density of the data and vastness of the gene set analyzed, this relatively low deviation
    # is sufficient to represent differential 3'UTR usage between cell clusters.
    cutoff = 50
    # Makes a hexbin plot using the average UTR deviation of each bin for coloration.
    plt.figure(figsize=(150, 150))
    plt.hexbin(x=tsne1, y=tsne2, C=utrs, gridsize=500, vmin=-cutoff, vmax=cutoff)
    plt.axis('off')
    tsne_figure_path = PAS_DATASET + "/figures/tsne_hex_utr_medians.png"
    plt.savefig(tsne_figure_path, bbox_inches='tight')
    plt.close()

    # Makes a hexbin density plot, to show the relative density contributing to the UTR deviations.
    plt.figure(figsize=(150, 150))
    plt.hexbin(x=tsne1, y=tsne2, gridsize=500, vmin=5, vmax=20, mincnt=1)
    plt.axis('off')
    kde_figure_path = PAS_DATASET + "/figures/tsne_hex_density.png"
    plt.savefig(kde_figure_path, bbox_inches='tight')
    plt.close()

    # Plots each cell age with a hexbin binned by average UTR deviation within that age.
    fig, axs = plt.subplots(1, 5, figsize=(50, 10), dpi=1000)
    for age_idx, age in enumerate(AGES):
        age_idxs = [idx for idx, entry in enumerate(ages) if entry == age]
        age_tsne1 = [tsne1[idx] for idx in age_idxs]
        age_tsne2 = [tsne2[idx] for idx in age_idxs]
        age_age_utr = [age_utrs[idx] for idx in age_idxs]
        axs[age_idx].hexbin(x=age_tsne1, y=age_tsne2, C=age_age_utr, gridsize=500, vmin=-cutoff, vmax=cutoff)
        axs[age_idx].set_title(age)
        axs[age_idx].axis('off')
    stsne_figure_path = PAS_DATASET + "/figures/age_tsne_hex_age_utr.png"
    plt.savefig(stsne_figure_path, bbox_inches='tight')
    plt.close()
    
    #  Plots each cell age with a hexbin binned by average UTR deviation within the full dataset.
    fig, axs = plt.subplots(1, 5, figsize=(50, 10), dpi=1000)
    for age_idx, age in enumerate(AGES):
        age_idxs = [idx for idx, entry in enumerate(ages) if entry == age]
        age_tsne1 = [tsne1[idx] for idx in age_idxs]
        age_tsne2 = [tsne2[idx] for idx in age_idxs]
        age_utr = [utrs[idx] for idx in age_idxs]
        axs[age_idx].hexbin(x=age_tsne1, y=age_tsne2, C=age_utr, gridsize=500, vmin=-cutoff, vmax=cutoff)
        axs[age_idx].set_title(age)
        axs[age_idx].axis('off')
    stsne_figure_path = PAS_DATASET + "/figures/age_tsne_hex_utr.png"
    plt.savefig(stsne_figure_path, bbox_inches='tight')
    plt.close()
    
    #  Plots each cell cluster with a hexbin binned by average UTR deviation within that cell cluster.
    fig, axs = plt.subplots(6, 7, figsize=(25, 25), dpi=1000)
    for cluster_num, cluster in enumerate(CELL_TYPES):
        cluster_row = cluster_num // 7
        cluster_col = cluster_num % 7
        cluster_idxs = [idx for idx, entry in enumerate(clusters) if entry == cluster]
        cluster_tsne1 = [stsne1[idx] for idx in cluster_idxs]
        cluster_tsne2 = [stsne2[idx] for idx in cluster_idxs]
        cluster_cluster_utr = [cluster_utrs[idx] for idx in cluster_idxs]
        axs[cluster_row, cluster_col].hexbin(x=cluster_tsne1, y=cluster_tsne2, C=cluster_cluster_utr,
                                             gridsize=300, vmin=-cutoff, vmax=cutoff)
        axs[cluster_row, cluster_col].set_title(cluster.replace('_', ' '), fontsize=10)
        axs[cluster_row, cluster_col].axis('off')
    # Blanks the remaining axes, as our dataset doesn't cover all 42 axes.
    for num in range(3, 7):
        axs[5, num].axis('off')
    stsne_figure_path = PAS_DATASET + "/figures/sub_tsne_hex_cluster_utr.png"
    plt.savefig(stsne_figure_path, bbox_inches='tight')
    plt.close()
    
    #  Plots each cell cluster with a hexbin binned by average UTR deviation within the full dataset.
    fig, axs = plt.subplots(6, 7, figsize=(25, 25), dpi=1000)
    for cluster_num, cluster in enumerate(CELL_TYPES):
        cluster_row = cluster_num // 7
        cluster_col = cluster_num % 7
        cluster_idxs = [idx for idx, entry in enumerate(clusters) if entry == cluster]
        cluster_tsne1 = [stsne1[idx] for idx in cluster_idxs]
        cluster_tsne2 = [stsne2[idx] for idx in cluster_idxs]
        cluster_utr = [utrs[idx] for idx in cluster_idxs]
        axs[cluster_row, cluster_col].hexbin(x=cluster_tsne1, y=cluster_tsne2, C=cluster_utr,
                                             gridsize=300, vmin=-cutoff, vmax=cutoff)
        axs[cluster_row, cluster_col].set_title(cluster.replace('_', ' '), fontsize=10)
        axs[cluster_row, cluster_col].axis('off')
    for num in range(3, 7):
        axs[5, num].axis('off')
    stsne_figure_path = PAS_DATASET + "/figures/sub_tsne_hex_utr.png"
    plt.savefig(stsne_figure_path, bbox_inches='tight')
    plt.close()


def bin_3d_distribution(x, y, z, weights, qual, bins, rounding=0, cutoff=1):
    """Turns three vectors of points into tessellated cubic bins, and then plots these bins as points.
    X, Y, and Z are the data vectors, with weights corresponding to data values.
    Qual is any qualitative array for each point. The function returns the most common qualitative value in each bin.
    Bins is a numeric of the number of bins in each dimension to use.
    Rounding is how to round the boundaries of the distribution.
    Cutoff is the number of entries a bin needs to not be discarded.
    Blows up data to 1000x for mathematical purposes, remember to alter your axes."""
    global_min = round(min(x + y + z), rounding) * 1000
    global_max = round(max(x + y + z), rounding) * 1000
    x = np.asarray(x) * 1000
    y = np.asarray(y) * 1000
    z = np.asarray(z) * 1000
    w = np.asarray(weights)
    a = np.asarray(qual)
    bin_boundary = (global_max - global_min) / bins
    binned_x = []
    binned_y = []
    binned_z = []
    binned_w = []
    binned_a = []
    x_count = 0
    x_min = global_min
    while len(w) > 0:
        x_count += 1
        print("X bin " + str(x_count) + " reached.")
        x_max = x_min + bin_boundary
        x_idxs = [idx for idx, coord in enumerate(x) if x_min <= coord < x_max]
        # Any time a range has no entries within it, the rest of the calculations are skipped.
        if len(x_idxs) < cutoff:
            x_min = x_max
            x = np.delete(x, x_idxs)
            y = np.delete(y, x_idxs)
            z = np.delete(z, x_idxs)
            w = np.delete(w, x_idxs)
            a = np.delete(a, x_idxs)
            continue
        # Subsets our array so we only search through entries that already satisfy our outer bins'
        # requirements. We do this for each of the arrays we search through, recursively.
        x_bin_y = y[[x_idxs]]
        x_bin_z = z[[x_idxs]]
        x_bin_w = w[[x_idxs]]
        x_bin_a = a[[x_idxs]]
        y_min = global_min
        while len(x_bin_w) > 0:
            y_max = y_min + bin_boundary
            y_idxs = [idx for idx, coord in enumerate(x_bin_y) if y_min <= coord < y_max]
            if len(y_idxs) < cutoff:
                y_min = y_max
                x_bin_y = np.delete(x_bin_y, y_idxs)
                x_bin_z = np.delete(x_bin_z, y_idxs)
                x_bin_w = np.delete(x_bin_w, y_idxs)
                x_bin_a = np.delete(x_bin_a, y_idxs)
                continue
            y_bin_z = x_bin_z[[y_idxs]]
            y_bin_w = x_bin_w[[y_idxs]]
            y_bin_a = x_bin_a[[y_idxs]]
            z_min = global_min
            while len(y_bin_w) > 0:
                z_max = z_min + bin_boundary
                z_idxs = [idx for idx, coord in enumerate(y_bin_z) if z_min <= coord < z_max]
                if len(z_idxs) < cutoff:
                    z_min = z_max
                    y_bin_z = np.delete(y_bin_z, z_idxs)
                    y_bin_w = np.delete(y_bin_w, z_idxs)
                    y_bin_a = np.delete(y_bin_a, z_idxs)
                    continue
                z_bin_w = y_bin_w[[z_idxs]]
                z_bin_a = y_bin_a[[z_idxs]]
                # Makes an X, Y, Z-centered bin and takes the average of the bin's weight values as the bin's
                # new 'weight' value.
                binned_x.append((x_max + x_min) / 2)
                binned_y.append((y_max + y_min) / 2)
                binned_z.append((z_max + z_min) / 2)
                binned_w.append(np.mean(z_bin_w))
                binned_a.append(sps.mode(z_bin_a)[0][0])
                # Everything below this clears binned entries from memory, saving time.
                z_min = z_max
                y_bin_z = np.delete(y_bin_z, z_idxs)
                y_bin_w = np.delete(y_bin_w, z_idxs)
            y_min = y_max
            x_bin_y = np.delete(x_bin_y, y_idxs)
            x_bin_z = np.delete(x_bin_z, y_idxs)
            x_bin_w = np.delete(x_bin_w, y_idxs)
        x_min = x_max
        x = np.delete(x, x_idxs)
        y = np.delete(y, x_idxs)
        z = np.delete(z, x_idxs)
        w = np.delete(w, x_idxs)
    return binned_x, binned_y, binned_z, binned_w, binned_a


def format_umap_data(refined_trajectories, refined_umap1, refined_umap2, refined_umap3, refined_umap_utr,
                     refined_umap_ages, formatted_trajectory_path):
    """Bins our UMAP data into cubic bins, and then writes the binned distribution to a text file for plotting in R."""
    with open(formatted_trajectory_path, 'wt') as trajectory_umap_out:
        # We choose these trajectories as they are the most differentiated trajectories in the global dataset.
        # This analysis can be done with any number of trajectories by simply adding their names to this list.
        for t_idx, traj in enumerate(['Neural_crest_melanocytes_trajectory', 'Neural_tube_and_notochord_trajectory',
                                      'Haematopoiesis_trajectory']):
            traj_idxs = [idx for idx, ref_traj in enumerate(refined_trajectories) if ref_traj == traj]
            traj_ref_umap1 = [refined_umap1[idx] for idx in traj_idxs]
            traj_ref_umap2 = [refined_umap2[idx] for idx in traj_idxs]
            traj_ref_umap3 = [refined_umap3[idx] for idx in traj_idxs]
            traj_ref_utr = [refined_umap_utr[idx] for idx in traj_idxs]
            traj_ref_age = [refined_umap_ages[idx] for idx in traj_idxs]
            binned_x, binned_y, binned_z, binned_w, binned_a = bin_3d_distribution(traj_ref_umap1, traj_ref_umap2,
                                                                                   traj_ref_umap3, traj_ref_utr,
                                                                                   traj_ref_age,
                                                                                   150, rounding=2, cutoff=2)
            # Writes our data into the output file.
            for bin_idx, weight in enumerate(binned_w):
                bin_x = str(binned_x[bin_idx])
                bin_y = str(binned_y[bin_idx])
                bin_z = str(binned_z[bin_idx])
                # This is already a string.
                bin_a = binned_a[bin_idx]
                data_for_r = '\t'.join([bin_x, bin_y, bin_z, str(weight), traj, bin_a]) + '\n'
                trajectory_umap_out.write(data_for_r)


def process_single_cell_data(formatted_trajectory_path):
    """Reads and formats our single-cell data for plotting."""
    # Reads in the single-cell data from our TSNE files generated by tsv_matrixlengthandisoformanalysis.py
    tsv_count = 0
    cell_data = []
    for tsv_file_path in glob2.glob(PAS_DATASET + '/tsv3/*.tsv'):
        with open(tsv_file_path, 'rt') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            for row in reader:
                if row[0] != 'NA':
                    cell_data.append(row)
        tsv_count += 1
        print("Cell TSV file " + str(tsv_count) + " processed.")

    # Generates arrays of data relevant to our TSNE plotting.
    clusters = [data[0] for data in cell_data]
    tsne1 = [float(data[1]) for data in cell_data]
    tsne2 = [float(data[2]) for data in cell_data]
    stsne1 = [float(data[4]) for data in cell_data]
    stsne2 = [float(data[5]) for data in cell_data]
    ages = [data[14] for data in cell_data]
    utrs = [float(data[15]) for data in cell_data]
    cluster_utrs = [float(data[16]) for data in cell_data]
    age_utrs = [float(data[19]) for data in cell_data]

    # Generates arrays of data relevant to our UMAP plotting.
    refined_umap1_raw = [data[11] for data in cell_data]
    refined_umap_idx = [idx for idx, umap in enumerate(refined_umap1_raw) if umap != 'NA']
    refined_umap_data = [cell_data[idx] for idx in refined_umap_idx]
    refined_trajectories = [data[6] for data in refined_umap_data]
    refined_umap1 = [float(data[11]) for data in refined_umap_data]
    refined_umap2 = [float(data[12]) for data in refined_umap_data]
    refined_umap3 = [float(data[13]) for data in refined_umap_data]
    refined_umap_utr = [float(data[17]) for data in refined_umap_data]
    refined_umap_ages = [data[14] for data in cell_data]

    generate_tsne_plots(tsne1, tsne2, utrs, ages, age_utrs, clusters, stsne1, stsne2, cluster_utrs)

    format_umap_data(refined_trajectories, refined_umap1, refined_umap2, refined_umap3, refined_umap_utr,
                     refined_umap_ages, formatted_trajectory_path)


def concatenate_matrices():
    """Takes the centered sparse matrices created by our 'centered' subscript and concatenates them, then saves them."""

    print("Concatenating processed data...")

    cells = []
    genes = []
    utrs = []
    counts = []

    # Simple concatenation of all of our fully processed data.
    for matrix_file in glob2.glob(PAS_DATASET + "/centered/*.h5"):
        with h5.File(matrix_file, "r") as input_file:
            cell = list(input_file["cells"])
            gene = list(input_file["genes"])
            utr = list(input_file["utrs"])
            count = list(input_file["counts"])
        cells = cells + cell
        genes = genes + gene
        utrs = utrs + utr
        counts = counts + count

    # We keep data sorted by cell index by default, since cell grouping is the most common way to sort the dataset.
    sorted_by_cell = sorted(zip(cells, genes, utrs, counts))
    cells_by_cell = [x[0] for x in sorted_by_cell]
    genes_by_cell = [x[1] for x in sorted_by_cell]
    utrs_by_cell = [x[2] for x in sorted_by_cell]
    counts_by_cell = [x[3] for x in sorted_by_cell]

    # Each of these objects are always added to the existing "overmatrix.h5" file, but the name of each data set is
    # different based on the input PAS dataset. This allows us to easily reference any of our created datasets.
    with h5.File("overmatrix.h5", "a") as out_file:
        out_file.create_dataset(PAS_DATASET + "_cells", data=cells_by_cell)
        out_file.create_dataset(PAS_DATASET + "_genes", data=genes_by_cell)
        out_file.create_dataset(PAS_DATASET + "_utrs", data=utrs_by_cell)
        out_file.create_dataset(PAS_DATASET + "_counts", data=counts_by_cell)

    print("Processed data concatenated.")


def main():
    # This will be passed our current directory from the bash pipeline script.
    os.chdir(sys.argv[2])
    sys.setrecursionlimit(10000)

    trans_to_gene_path = "trans_to_gene.txt"
    if os.path.exists(trans_to_gene_path):
        print("Transcript to gene association already complete.")
    else:
        print("Choosing most expressed transcripts...")
        choose_transcripts(trans_to_gene_path)

    # These variables are all referenced and/or modified by multiple functions, so they are set as globals.
    global GENE_NAMES
    global APPROVED_TRANSCRIPTS
    with open(TRANSCRIPT_FILE_PATH, 'rt') as transcript_file:
        APPROVED_TRANSCRIPTS = [line.rstrip('\n') for line in transcript_file]

    cell_data_path = "cell_data_dict.pkl"
    names_info_path = "names_by_id.pkl"
    gene_data_path = "gene_data_dict.pkl"
    if os.path.exists(cell_data_path):
        print("Cell data found.")
    else:
        print("Cell data file missing.")
        generate_cell_data(cell_data_path)
    if os.path.exists(names_info_path):
        print("Data already initiated")
    else:
        print("Initiating data files...")
        initiate_data_files(names_info_path)
    if os.path.exists(gene_data_path):
        print("Gene data dict found")
    else:
        print("Gene data dict missing.")
        generate_gene_data_dict(gene_data_path)

    # If we enter a PAS isoform dataset to process, this will specify the folder that it's contained within.
    # In the case that we don't provide a dataset, the default dataset is "raw" which isn't actually a dataset at all.
    global PAS_DATASET
    if PAS_DATASET != "raw":
        pas_venn_path = PAS_DATASET + "/figures/pas_dataset_overlap.pdf"
        if os.path.exists(pas_venn_path):
            print("PAS information found.")
        else:
            print("Generating PAS information...")
            generate_pas_files_2()
            print("Analyzing PAS database overlap...")
            generate_pas_venn(pas_venn_path)

    gene_average_path = PAS_DATASET + "/gene_average_dict.pkl"
    if os.path.exists(gene_average_path):
        print("Gene average dict found.")
    else:
        print("Gene average dict missing.")
        print("Running raw data analysis...")
        # This builds our initial data matrices so we can calculate gene averages.
        coverage_path = PAS_DATASET + '/coverage/'
        if not os.path.exists(coverage_path):
            os.mkdir(coverage_path)
        parallel_job_submission("raw")
        print("Analyzing gene coverage...")
        analyze_gene_coverage()
        print("Calculating average gene utr lengths...")
        calculate_averages(gene_average_path)

    pas_matrix_path = PAS_DATASET + '/pas_matrices.h5'
    #if os.path.exists(pas_matrix_path):
    #    print("PAS expression data found.")
    #else:
    #    print("Analyzing pas usage...")
    analyze_pas_usage()

    transcript_count_path = PAS_DATASET + "/transcript_counts_pasoverlaps.txt"
    if os.path.exists(transcript_count_path):
        print("Transcript counts found.")
    else:
        print("Generating transcript counts...")
        parallel_job_submission("count")
        count_transcripts(transcript_count_path)
    gene_expression_path = PAS_DATASET + "/gene_expr_data.txt"
    if os.path.exists(gene_expression_path):
        print("Gene expression data formatted")
    else:
        print("Formatting gene expression data...")
        format_expression_data(gene_expression_path)

    print("Centering data...")
    tsv_path = PAS_DATASET + "/tsv/"
    heatmap_array_path = PAS_DATASET + "/heatmap/"
    pas_array_path = PAS_DATASET + "/pasarrays/"
    heatmap_path = PAS_DATASET + "/figures/gene_age_heatmap.pdf"
    gene_stats_path = PAS_DATASET + "/gene_stats.txt"
    pas_stats_path = PAS_DATASET + "/pas_stats_cluster.txt"
    formatted_trajectory_path = PAS_DATASET + "/binned_umap_data.txt"
    if not os.path.exists(tsv_path):
        os.mkdir(tsv_path)
    if not os.path.exists(heatmap_array_path):
        os.mkdir(heatmap_array_path)
    if not os.path.exists(pas_array_path):
        os.mkdir(pas_array_path)
    # This will not only center our data but also make our cell tsvs. We can concatenate them as we see fit.
    parallel_job_submission("centered")
    if os.path.exists(heatmap_path):
        print("Heatmaps already generated.")
    else:
        print("Generating heatmaps...")
        generate_heatmaps(heatmap_path)
    if os.path.exists(gene_stats_path):
        print("Gene stats already generated.")
    else:
        print("Calculating gene stats...")
        calculate_transcript_stats(gene_stats_path)
    if os.path.exists(pas_stats_path):
        print("PAS stats already generated.")
    else:
        print("Calculating pas stats...")
        calculate_pas_stats(pas_stats_path)
    # There's no followup concatenation function to this, as we'd prefer to just concatenate these tsvs in bash.
    parallel_job_submission("tsv")
    if os.path.exists(formatted_trajectory_path):
        print("Single cell plots already generated.")
    else:
        print("Formatting single cell data...")
        process_single_cell_data(formatted_trajectory_path)
    print("All data generated.")

    # This is only to concatenate all of the generated data into a cell by gene numpy matrix stored in an h5 file.
    # This process takes several days on my cluster, and may take longer on yours. Run at your own leisure.
    # with h5.File("overmatrix.h5", "a") as overmatrix_in:
    #     if PAS_DATASET + "_counts" in list(overmatrix_in.keys()):
    #         print("Concatenated data already stored.")
    #     else:
    #         concatenate_matrices()


if __name__ == "__main__":
    main()

