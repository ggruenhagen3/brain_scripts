import scipy.sparse as sparse
import numpy as np
import pandas

data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_data_mat.npz")
gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_rownames.csv").iloc[:,
              1].to_numpy()
cond_labels = pandas.read_csv(
    "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_real_cond_labels.csv").iloc[:,
              1].to_numpy()

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
    # Find Node Strength
    ns = cor.sum(axis=1)
    print(id(data_mat))
    return ns

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
    print("Loaded Test")

if __name__ == '__main__':
    main()