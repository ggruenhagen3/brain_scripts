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
import grn_test

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
    # mat = data_mat[:, this_idx]
    # Find Correlations
    print("Finding correlations")
    cor = pandas.DataFrame(data = sparse_corrcoef(data_mat[:, this_idx].todense()), index = gene_labels, columns = gene_labels)
    print("Done w/ corr")
    # Find Node Strength
    ns = cor.sum(axis=1)
    print(id(data_mat))
    return ns


def sparse_corrcoef(A, B=None):
    print("In corr calc")
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
    print("Before corr error")
    coeffs = C / np.sqrt(np.outer(d, d))
    print("After corr error")
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
    #Globals
    global data_mat
    global gene_labels
    global cond_labels
    data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_data_mat.npz")
    gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_rownames.csv").iloc[:,
                  1].to_numpy()
    cond_labels = pandas.read_csv(
        "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_real_cond_labels.csv").iloc[:,
                  1].to_numpy()
    # Read Inputs
    perm_num, num_perm = parseArgs()
    # Read BB data
    # data_mat = pickle.load(open("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_data_mat.pickle", "rb")) # data matrix
    # Set random seed so all the permutations are different
    print("Seed = " + str(perm_num))
    random.seed(perm_num)
    # Permute BHVE and CTRL labels
    print("Permuting Data " + str(num_perm) + " times.")
    mat_idx = permuteLabels(num_perm)
    print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")
    # Create BHVE and CTRL Matrices, Find Correlations, Find Node Strengths and Find NodeStrength Differences
    # this_result = pool.map(corAndNodeStrength, [perm_label for perm_label in perm_labels])
    mat_idx3 = {'B0': [0, 1, 2, 3, 5], 'C0': [4, 6, 10, 11, 15], 'B1': [1, 3, 4, 7, 8]}
    print("Created Mat 3")
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        ns_dict = pool.map(corAndNodeStrength, mat_idx3.values())
    # df=pd.DataFrame.from_dict(d,orient='index').transpose()
    # test = corAndNodeStrength(cond_labels)

# def testFun(values):
#     print(values)
#
# mat_idx2 = {}
# for key in mat_idx.keys():
#     mat_idx2[key] = mat_idx[key][0:5]


if __name__ == '__main__':
    main()
