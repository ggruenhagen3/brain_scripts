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

def parseArgs():
    parser = argparse.ArgumentParser(description='Shuffle BHVE and CTRL, find correlations for both, node strength for both and node strength difference.')
    parser.add_argument('perm_num', metavar='perm_num', type = int, help='The current permutation number. This is used for the seed.')
    parser.add_argument('num_perm', metavar='num_perm', type = int, help='The number of permutations to complete.')
    args = parser.parse_args()
    return args.perm_num, args.num_perm

# def corAndNodeStrength(perm_label, data_mat = data_mat, gene_labels = gene_labels):
def corAndNodeStrength(perm_label):
    """
    For 1 permutation, create BHVE and CTRL Matrices, find correlations, find node strengths and find node strength differences.
    :param perm_label: Labels of condition (BHVE and CTRL) for this permutation
    :return bhve_minus_ctrl_ns: bhve node strength - control node strength
    """
    # Find BHVE and CTRL cells
    bhve_idx = []
    ctrl_idx = []
    for i in range(0, len(perm_label)):
        this_label = perm_label[i]
        if this_label == "BHVE":
            bhve_idx.append(i)
        else:
            ctrl_idx.append(i)
    # Create BHVE and CTRL Matrices
    bhve_mat = data_mat[:, bhve_idx]
    ctrl_mat = data_mat[:, ctrl_idx]
    # Find Correlations
    print("Finding BHVE correlations")
    bhve_cor = pandas.DataFrame(data = sparse_corrcoef(bhve_mat.todense()), index = gene_labels, columns = gene_labels)
    print("Finding CTRL correlations")
    ctrl_cor = pandas.DataFrame(data = sparse_corrcoef(ctrl_mat.todense()), index = gene_labels, columns = gene_labels)
    # Find Node Strength
    bhve_ns = bhve_cor.sum(axis=1)
    ctrl_ns = ctrl_cor.sum(axis=1)
    bhve_minus_ctrl_ns = bhve_ns - ctrl_ns
    print(bhve_ns[0:5])
    print(ctrl_ns[0:5])
    print(bhve_minus_ctrl_ns[0:5])
    return(bhve_minus_ctrl_ns)


def sparse_corrcoef(A, B=None):
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

def main():
    # Start the timer
    start_time = time.perf_counter()
    # Read Inputs
    perm_num, num_perm = parseArgs()
    # Read BB data
    # data_mat = pickle.load(open("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_data_mat.pickle", "rb")) # data matrix
    data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_data_mat.npz")
    gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_rownames.csv").iloc[:,1].to_numpy()
    cond_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_real_cond_labels.csv").iloc[:,1].to_numpy()
    # Set random seed so all the permutations are different
    print("Seed = " + str(perm_num))
    random.seed(perm_num)
    # Permute BHVE and CTRL labels
    print("Permuting Data" + str(num_perm) + " times.")
    perm_labels = []
    for i in range(0, num_perm):
        perm_labels.append(random.shuffle(cond_labels))
    print(f"Done Permuting. Current Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")
    # Create BHVE and CTRL Matrices, Find Correlations, Find Node Strengths and Find NodeStrength Differences
    # pool = multiprocessing.Pool(num_perm)
    # this_result = pool.map(corAndNodeStrength, [perm_label for perm_label in perm_labels])
    test = corAndNodeStrength(cond_labels)


if __name__ == '__main__':
    main()
