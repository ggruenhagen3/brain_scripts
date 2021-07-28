# Load packages
import scipy.sparse as sparse
import pandas
import numpy as np
import time
import random
import argparse
import multiprocessing
from pathlib import Path
from itertools import repeat

global data_mat
global gene_labels
global cluster_labels
global do_abs
global cluster_set

def parseArgs():
    parser = argparse.ArgumentParser(description='Shuffle cluster labels and see if a gene has significantly greater node strength in the real vs perm.')
    parser.add_argument("dataset", metavar="dataset", type=str, help="Dataset (ct, cj, mi, mie, mim, hm)")
    parser.add_argument('perm_num', metavar='perm_num', type = int, help='The current permutation number. This is used for the seed.')
    parser.add_argument('num_perm', metavar='num_perm', type = int, help='The number of permutations to complete.')
    parser.add_argument("-o", "--output_folder", help="Output Folder", nargs="?",
                        default="/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/py_ns/",
                        const="/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/py_ns/")
    parser.add_argument("-n", "--no_perm", help="Do no permutations?", action="store_true")
    parser.add_argument("-a", "--do_abs", help="Take the absolute value of the correlations?", action="store_true")
    args = parser.parse_args()
    return args.perm_num, args.num_perm, args.output_folder, args.no_perm, args.dataset, args.do_abs

def myShuffle(this_list):
    """
    Shuffle and return the list. (I made this function bc I don't like the how regular shuffle returns the list).
    :param this_list: input list (unshuffled)
    :return this_list: shuffled list
    """
    random.shuffle(this_list)
    return(this_list)

def corAndNodeStrength(this_idx):
    """
    Given idexes of cells, create a matrix, find correlations, and find node strengths
    :param i: i
    :param cluster: Cluster
    :return ns: Node Strength
    """
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

def permuteLabels(num_perm):
    """
    Permute cluster labels and find the indexes of the matrix belong to each.
    :param num_perm: Number of permutations
    :return mat_idx: dictionary of lists of column indexes that belong to BHVE and CTRL for each permutation
    """
    label_dict = {} # key is permutation #, value is a smaller dictionary
    for i in range(0, num_perm):
        perm_label = myShuffle(cluster_labels)
        small_dict = {}  # key is the cluster and value is indexes of the cluster
        for cluster in cluster_set:
            # Idx of labels equal to this cluster
            small_dict[cluster] = np.hstack(np.argwhere(perm_label == cluster))
        label_dict[i] = small_dict
    return label_dict

def main():
    # Start the timer
    start_time = time.perf_counter()

    # Read Inputs
    global data_mat
    global gene_labels
    global cluster_labels
    global do_abs
    global cluster_set
    perm_num, num_perm, output_folder, no_perm, dataset, do_abs = parseArgs()

    # Check Dataset Input
    if dataset == "ct":
        data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/tj_data_mat.npz")
        gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/tj_names.csv").iloc[:, 1].to_numpy()
        cluster_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/tj_clusters.csv").iloc[:, 1].to_numpy()
    elif dataset == "cj":
        data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/jaw_data_mat.npz")
        gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/jaw_names.csv").iloc[:, 1].to_numpy()
        cluster_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/jaw_clusters.csv").iloc[:, 1].to_numpy()
    elif dataset == "mi":
        data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/incsr_data_mat.npz")
        gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/incsr_names.csv").iloc[:, 1].to_numpy()
        cluster_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/incsr_clusters.csv").iloc[:, 1].to_numpy()
    elif dataset == "mie":
        data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/mie_data_mat.npz")
        gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/mie_names.csv").iloc[:, 1].to_numpy()
        cluster_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/mie_clusters.csv").iloc[:, 1].to_numpy()
    elif dataset == "mim":
        data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/im_data_mat.npz")
        gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/im_names.csv").iloc[:, 1].to_numpy()
        cluster_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/im_clusters.csv").iloc[:, 1].to_numpy()
    elif dataset == "hm":
        data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/hm_data_mat.npz")
        gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/hm_names.csv").iloc[:, 1].to_numpy()
        cluster_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/hm_clusters.csv").iloc[:, 1].to_numpy()
    else:
        print("Invalid dataset, please select one of ct, cj, mi, mie, mim, hm.")
        return

    cluster_set = list(set(cluster_labels))

    # Change folder name based on input
    base_name = dataset + "_perm_" + str(perm_num)
    if no_perm:
        base_name = dataset + "_real"
    if do_abs:
        base_name = base_name + "_abs"

    # Set random seed so all the permutations are different
    print("Seed = " + str(perm_num))
    random.seed(perm_num)

    # Permute Cluster Labels
    if no_perm:
        print("Not permuting Data")
        mat_idx = {}
        small_dict = {}  # key is the cluster and value is indexes of the cluster
        for cluster in cluster_set:
            # Idx of labels equal to this cluster
            small_dict[cluster] = np.hstack(np.argwhere(cluster_labels == cluster))
        mat_idx[0] = small_dict
    else:
        print("Permuting Data " + str(num_perm) + " times.")
        mat_idx = permuteLabels(num_perm)
        print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")

    # Set up dataframes that store results
    all_cluster_df = {}
    for cluster in cluster_set:
        all_cluster_df[cluster] = pandas.DataFrame(index = gene_labels, columns = list(range(1,num_perm+1)))

    # Find Correlations and Node Strengths
    print("Finding Correlations")
    for i in range(0, num_perm):
        print("Perm: " + str(i))
        perm_start = time.perf_counter()
        with multiprocessing.Pool(len(cluster_set)) as pool:
            this_idx_list = mat_idx[i].values()
            pool_ns = pool.map(corAndNodeStrength, this_idx_list)
            # pool_ns = pool.starmap(corAndNodeStrength, zip(repeat(i), cluster_set))
            # for j in range(0, len(cluster_set)):
            #     all_cluster_df[cluster_set[j]][i+1] = pool_ns[j]
        print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - perm_start:0.4f} seconds")

    # Write results for each cluster to a file
    for cluster in cluster_set:
        Path(output_folder + "/" + base_name + "/").mkdir(parents=True, exist_ok=True)
        all_cluster_df[cluster].to_csv(output_folder + "/" + base_name + "/" + cluster + "_" + str(perm_num) + ".csv")

if __name__ == '__main__':
    main()