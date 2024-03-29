# Load packages
import pickle
import scipy
import scipy.sparse as sparse
import pandas
import numpy as np
import time
import random
import argparse
import h5py
import multiprocessing
from itertools import repeat

global data_mat
global gene_labels
global cond_labels

def parseArgs():
    parser = argparse.ArgumentParser(description='Shuffle BHVE and CTRL, find correlations for both, node strength for both and node strength difference.')
    parser.add_argument('perm_num', metavar='perm_num', type = int, help='The current permutation number. This is used for the seed.')
    parser.add_argument('num_perm', metavar='num_perm', type = int, help='The number of permutations to complete.')
    parser.add_argument("-c", "--cluster15", help="15 Cluster to subset data by", nargs="?", type=int, default=-1, const=-1)
    parser.add_argument("-b", "--cluster53", help="53 Cluster to subset data by", nargs="?", type=int, default=-1, const=-1)
    parser.add_argument("-g", "--gene", help="Find correlations only in gene positive cells", nargs="?", type=str, default="", const="")
    parser.add_argument("-o", "--output_folder", help="Output Folder", nargs="?",
                        default="/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/py_ns/",
                        const="/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/py_ns/")
    parser.add_argument("-n", "--no_perm", help="Do no permutations?", action="store_true")
    parser.add_argument("-r", "--cor_only", help="Only find correlations in the data?", action="store_true")
    parser.add_argument("-a", "--do_abs", help="Take the absolute value of the correlations?", action="store_true")
    parser.add_argument("-s", "--sum_ns", help="Sum absolute value of node strength difference?", action="store_true")
    parser.add_argument("-m", "--replicate_match", help="Does the NS Dif match for all replicates?", action="store_true")
    parser.add_argument("-i", "--ieg", help="Run on cells that express at least X IEGs.", nargs="?", type=int, default=0, const=0)
    parser.add_argument("-e", "--no_bvc", help="Find cor in all cells instead of BHVE vs CTRL?", action="store_true")
    args = parser.parse_args()
    return args.perm_num, args.num_perm, args.cluster15, args.cluster53, args.gene, args.output_folder, args.no_perm, args.cor_only, args.do_abs, args.sum_ns, args.replicate_match, args.ieg, args.no_bvc


def corOnlyAndWrite(this_idx, output_path):
    """
    Given idexes of cells, create a matrix and find correlations only
    :param this_idx: Indexes of columns
    :param output_path: Output path of h5 correlation matrix file
    :return success: Function completed? True/False
    """
    cor = pandas.DataFrame(data=sparse_corrcoef(data_mat[:, this_idx].todense()), index=gene_labels, columns=gene_labels)
    if do_abs:
        print("Taking absolute value of correlations")
        cor = cor.abs()
    else:
        print("NOT taking absolute value of correlations. Using raw values.")
    h5f = h5py.File(output_path, 'w')
    h5f.create_dataset('name', data=cor)
    h5f.close()
    return True

def corAndNodeStrength(this_idx):
    """
    Given idexes of cells, creat a matrix, find correlations, and find node strengths
    :param this_idx: Indexes of columns
    :return ns: Node Strength
    """
    # Create BHVE and CTRL Matrices
    # Find Correlations
    cor = pandas.DataFrame(data = sparse_corrcoef(data_mat[:, this_idx].todense()), index = gene_labels, columns = gene_labels)
    if do_abs:
        print("Taking absolute value of correlations")
        cor = cor.abs()
    else:
        print("NOT taking absolute value of correlations. Using raw values.")
    # Find Node Strength
    ns = cor.sum(axis=1)
    return ns


def sparse_corrcoef(A, B=None):
    """
    Find correlations in sparse matrix
    """
    if B is not None:
        A = sparse.vstack((A, B), format='csr')
    A = A.astype(np.float64)
    n = A.shape[1]
    # Compute the covariance matrix
    rowsum = A.sum(1)
    centering = rowsum.dot(rowsum.T.conjugate()) / n
    C = (A.dot(A.T.conjugate()) - centering) / (n - 1)
    # The correlation coefficients are given by
    # C_{i,j} / sqrt(C_{i} * C_{j})
    d = np.diag(C)
    coeffs = C / np.sqrt(np.outer(d, d))
    return coeffs

def myShuffle(this_list):
    """
    Shuffle and return the list. (I made this function bc I don't like the how regular shuffle returns the list).
    :param this_list: input list (unshuffled)
    :return this_list: shuffled list
    """
    random.shuffle(this_list)
    return(this_list)

def permuteLabels(num_perm):
    """
    Permute BHVE and CTRL condition labels and find the indexes of the matrix belong to each.
    :param num_perm: Number of permutations
    :return mat_idx: dictionary of lists of column indexes that belong to BHVE and CTRL for each permutation
    """
    mat_idx = {} # key is name of cond + perm, value is columns that belong to that simulated sample
    for i in range(0, num_perm):
        perm_label = myShuffle(cond_labels)
        bhve_idx = []
        ctrl_idx = []
        for j in range(0, len(perm_label)):
            this_label = perm_label[j]
            if this_label == "BHVE":
                bhve_idx.append(j)
            else:
                ctrl_idx.append(j)
        mat_idx["B" + str(i)] = bhve_idx
        mat_idx["C" + str(i)] = ctrl_idx
    return mat_idx


def main():
    # Start the timer
    start_time = time.perf_counter()

    # Read Inputs
    global do_abs
    perm_num, num_perm, cluster15, cluster53, gene, output_folder, no_perm, cor_only, do_abs, sum_ns, replicate_match, ieg, no_bvc = parseArgs()

    # Read BB data
    global data_mat
    global gene_labels
    global cond_labels
    data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_data_mat.npz")
    gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_rownames.csv").iloc[:, 1].to_numpy()
    cond_labels = pandas.read_csv(
        "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_real_cond_labels.csv").iloc[:, 1].to_numpy()
    cluster15_labels = pandas.read_csv(
        "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_seuratclusters15.csv").iloc[:, 0].to_numpy()
    cluster53_labels = pandas.read_csv(
        "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_seuratclusters53.csv").iloc[:, 0].to_numpy()
    bb_metadata = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_metadata.csv")
    ieg_score = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_ieg_score.csv").iloc[:,1].to_numpy()

    # Subset by gene if necessary
    if gene != "":
        print("Subsetting by gene: " + gene)
        gene_idx = np.where(gene_labels == gene)[0][0]  # row number of the gene
        gene_pos_idx = np.nonzero(data_mat[gene_idx])[1] # gene positive cells
        data_mat = data_mat[:, gene_pos_idx]
        cond_labels = cond_labels[gene_pos_idx]

    # Subset by IEG Score if necessary
    if ieg > 0:
        ieg_idx = np.where(ieg_score >= ieg)[0]
        data_mat = data_mat[:, ieg_idx]
        cond_labels = cond_labels[ieg_idx]

    # Change file name based on input
    base_name = "perm_"
    if no_perm:
        base_name = "real_"
    if do_abs:
        base_name = base_name + "abs_"
    if gene != "":
        base_name = base_name + str(gene)
    if replicate_match:
        base_name = base_name + "replicate"
    if no_bvc:
        base_name = base_name + "no_bvc_"
    if ieg > 0:
        base_name = base_name + "ieg" + str(ieg)

    # Subset by cluster if necessary
    if cluster15 != -1:
        print("Subsetting on 15 cluster level for cluster " + str(cluster15))
        base_name = base_name + "cluster15_" + str(cluster15)
        data_mat = data_mat[:, np.flatnonzero(cluster15_labels == cluster15)]
        cond_labels = cond_labels[np.flatnonzero(cluster15_labels == cluster15)]
    else:
        if cluster53 != -1:
            print("Subsetting on 53 cluster level for cluster " + str(cluster53))
            base_name = base_name + "cluster53_" + str(cluster53)
            data_mat = data_mat[:, np.flatnonzero(cluster53_labels == cluster53)]
            cond_labels = cond_labels[np.flatnonzero(cluster53_labels == cluster53)]
        else:
            print("Not subsetting by any clusters.")

    # If necessary, find the correlations only (and don't do node strength stuff) and do no permutations
    if cor_only:
        print("Finding Correlations Only")
        base_name = base_name + "_cor"
        bhve_output_path = output_folder + "/" + base_name + "_bhve.h5"
        ctrl_output_path = output_folder + "/" + base_name + "_ctrl.h5"
        no_bvc_output_path = output_folder + "/" + base_name + "_all.h5"
        if no_bvc:
            print("Find Cor in ALL Cells, Not BHVE vs CTRL.")
            all_cor_success = corOnlyAndWrite(range(0, data_mat.shape[1]), no_bvc_output_path)
        else:
            bhve_cor_success = corOnlyAndWrite(np.flatnonzero(cond_labels == "BHVE"), bhve_output_path)
            ctrl_cor_success = corOnlyAndWrite(np.flatnonzero(cond_labels == "CTRL"), ctrl_output_path)
        print("Done")
        return

    # Set random seed so all the permutations are different
    print("Seed = " + str(perm_num))
    random.seed(perm_num)

    # Permute BHVE and CTRL labels
    if no_perm:
        print("Not permuting Data")
        mat_idx = {}
        mat_idx["B" + str(0)] = np.flatnonzero(cond_labels == "BHVE")
        mat_idx["C" + str(0)] = np.flatnonzero(cond_labels == "CTRL")
    else:
        print("Permuting Data " + str(num_perm) + " times.")
        mat_idx = permuteLabels(num_perm)
        print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")

    # Subset by cells if necessary
    if replicate_match:
        print("Checking to see all replicates have a matching direction for NS Dif")
        num_perm = 5
        mat_idx = {}
        mat_idx["B" + str(0)] = np.flatnonzero(bb_metadata[["sample"]] == "b1")
        mat_idx["C" + str(0)] = np.flatnonzero(bb_metadata[["sample"]] == "c1")
        mat_idx["B" + str(1)] = np.flatnonzero(bb_metadata[["sample"]] == "b2")
        mat_idx["C" + str(1)] = np.flatnonzero(bb_metadata[["sample"]] == "c2")
        mat_idx["B" + str(2)] = np.flatnonzero(bb_metadata[["sample"]] == "b3")
        mat_idx["C" + str(2)] = np.flatnonzero(bb_metadata[["sample"]] == "c3")
        mat_idx["B" + str(3)] = np.flatnonzero(bb_metadata[["sample"]] == "b4")
        mat_idx["C" + str(3)] = np.flatnonzero(bb_metadata[["sample"]] == "c4")
        mat_idx["B" + str(4)] = np.flatnonzero(bb_metadata[["sample"]] == "b5")
        mat_idx["C" + str(4)] = np.flatnonzero(bb_metadata[["sample"]] == "c5")

    # Create BHVE and CTRL Matrices, Find Correlations, Find Node Strengths and Find NodeStrength Differences
    # mat_idx3 = {'B0': [0, 1, 2, 3, 5], 'C0': [4, 6, 10, 11, 15], 'B1': [1, 3, 4, 7, 8]}
    print("Finding Correlations")
    ns_dict = {}
    for i in range(0, num_perm):
        print("Start Pair: " + str(i))
        with multiprocessing.Pool(2) as pool:
            pool_ns = pool.map(corAndNodeStrength, [mat_idx['B' + str(i)], mat_idx['C' + str(i)]])
            if no_perm and not replicate_match:
                ns_dict["B"] = pool_ns[0]
                ns_dict["C"] = pool_ns[1]
                ns_dict["Dif"] = pool_ns[0] - pool_ns[1]
                ns_dict["Dif_Abs"] = np.absolute(pool_ns[0] - pool_ns[1])
            else:
                ns_dict[i] = pool_ns[0] - pool_ns[1]
            # Find Genes that are in at least 5 cells in BHVE and CTRL
            num_cells_bhve = np.array(data_mat[:, mat_idx['B' + str(i)]].astype(bool).sum(axis=1))[:, 0]
            num_cells_ctrl = np.array(data_mat[:, mat_idx['C' + str(i)]].astype(bool).sum(axis=1))[:, 0]
            gene_idx_min_cells = np.where(np.logical_and(num_cells_bhve >= 5, num_cells_ctrl >= 5))[0]
            # Find Sum of NS and Sum of Absolute NS for those genes
            sum_ns_dif = np.sum(pool_ns[0][gene_idx_min_cells] - pool_ns[1][gene_idx_min_cells])
            sum_abs_ns_dif = np.sum(np.absolute(pool_ns[0][gene_idx_min_cells] - pool_ns[1][gene_idx_min_cells]))
            if sum_ns:
                ns_dict[i] = sum_abs_ns_dif
            print("Sums for genes with at least 5 cells in BHVE and CTRL")
            print("Sum of NS Dif " + str(sum_ns_dif) )
            print("Sum of Abs NS Dif " + str(sum_abs_ns_dif))
        print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")
    perm_ns_dif = pandas.DataFrame.from_dict(ns_dict,orient='index').transpose()
    if replicate_match:
        perm_ns_dif.columns = ["NS_Dif_1", "NS_Dif_2", "NS_Dif_3", "NS_Dif_4", "NS_Dif_5"]
    perm_ns_dif.to_csv(output_folder + "/" + base_name + "_" + str(num_perm) + ".csv")

if __name__ == '__main__':
    main()
