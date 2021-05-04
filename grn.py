# Load packages
import pickle
import scipy
import scipy.sparse as sparse
import pandas
import numpy as np
import time
import random
import argparse
import multiprocessing
from itertools import repeat

global data_mat
global gene_labels
global cond_labels

def parseArgs():
    parser = argparse.ArgumentParser(description='Shuffle BHVE and CTRL, find correlations for both, node strength for both and node strength difference.')
    parser.add_argument('perm_num', metavar='perm_num', type = int, help='The current permutation number. This is used for the seed.')
    parser.add_argument('num_perm', metavar='num_perm', type = int, help='The number of permutations to complete.')
    args = parser.parse_args()
    return args.perm_num, args.num_perm


def corAndNodeStrength(this_idx):
    """
    Given idexes of cells, creat a matrix, find correlations, and find node strengths
    :param this_idx: Indexes of columns
    :return ns: Node Strength
    """
    # Create BHVE and CTRL Matrices
    # Find Correlations
    cor = pandas.DataFrame(data = sparse_corrcoef(data_mat[:, this_idx].todense()), index = gene_labels, columns = gene_labels)
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
    perm_num, num_perm = parseArgs()
    # Read BB data
    global data_mat
    global gene_labels
    global cond_labels
    data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_data_mat.npz")
    gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_rownames.csv").iloc[:,
                  1].to_numpy()
    cond_labels = pandas.read_csv(
        "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_real_cond_labels.csv").iloc[:,
                  1].to_numpy()
    # Set random seed so all the permutations are different
    print("Seed = " + str(perm_num))
    random.seed(perm_num)
    # Permute BHVE and CTRL labels
    print("Permuting Data " + str(num_perm) + " times.")
    mat_idx = permuteLabels(num_perm)
    print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")
    # Create BHVE and CTRL Matrices, Find Correlations, Find Node Strengths and Find NodeStrength Differences
    # mat_idx3 = {'B0': [0, 1, 2, 3, 5], 'C0': [4, 6, 10, 11, 15], 'B1': [1, 3, 4, 7, 8]}
    print("Finding Correlations")
    ns_dict = {}
    for i in range(0, num_perm):
        print("Start Pair: " + str(i))
        with multiprocessing.Pool(2) as pool:
            pool_ns = pool.map(corAndNodeStrength, [mat_idx['B' + str(i)], mat_idx['C' + str(i)]])
            ns_dict[i] = pool_ns[0] - pool_ns[1]
        print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")
    perm_ns_dif = pandas.DataFrame.from_dict(ns_dict,orient='index').transpose()
    perm_ns_dif.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/py_ns" + str(perm_num))

if __name__ == '__main__':
    main()
