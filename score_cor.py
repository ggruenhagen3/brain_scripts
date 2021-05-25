# Load packages
import scipy.sparse as sparse
import pandas
import numpy as np
import random
import time
import argparse
import multiprocessing
from itertools import repeat

def parseArgs():
    parser = argparse.ArgumentParser(description='Find correlations and their differences in BHVE and CTRL for each gene with a score for a list of genes.')
    parser.add_argument('gene_list', metavar='gene_list', type = str, help='List of genes to create a score for. Options are ieg, neurogen, and prog.')
    parser.add_argument('num_perm', metavar='num_perm', type = int, help='The number of permutations to complete.')
    args = parser.parse_args()
    return args.gene_list, args.num_perm

def findCor(gene_idx, subset_idx = None):
    cor = 0
    if subset_idx is not None:
        cor = np.corrcoef(score[subset_idx], data_mat[gene_idx, subset_idx].todense())[0,1]
    else:
        cor = np.corrcoef(score, data_mat[gene_idx, ].todense())[0,1]
    return cor

def myShuffle(this_list):
    """
    Shuffle and return the list. (I made this function bc I don't like the how regular shuffle returns the list).
    :param this_list: input list (unshuffled)
    :return this_list: shuffled list
    """
    random.shuffle(this_list)
    return(this_list)

def singleRun():
    # Create Score
    mat = data_mat[:, :]
    nonzero_mat = mat.nonzero()
    mat[nonzero_mat[0], nonzero_mat[1]] = 1
    global score
    score = np.array(mat[score_genes_idx,].sum(axis=0)) / np.array(mat.sum(axis=0))
    score = score[0]

    # Find the correlated genes in bulk in both BHVE and CTRL
    bulk_res = {}
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        this_res = pool.starmap(findCor, zip(range(0, len(gene_labels)), repeat(bhve_idx, len(gene_labels))))
        bulk_res["BHVE"] = this_res
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        this_res = pool.starmap(findCor, zip(range(0, len(gene_labels)), repeat(ctrl_idx, len(gene_labels))))
        bulk_res["CTRL"] = this_res

    # Save the bulk results to a dataframe
    bulk_df = pandas.DataFrame(bulk_res, index=gene_labels)
    bulk_df['Dif'] = bulk_df['BHVE'] - bulk_df['CTRL']

    # 15 cluster level
    clust15_res = {}
    clust15_dif_res = {}
    for clust15 in range(0, 15):
        print(clust15)
        clust_idx = np.where(cluster15_labels == clust15)[0]
        clust_bhve_idx = pandas.Series(bhve_idx).isin(clust_idx)
        clust_bhve_idx = bhve_idx[clust_bhve_idx]
        clust_ctrl_idx = pandas.Series(ctrl_idx).isin(clust_idx)
        clust_ctrl_idx = ctrl_idx[clust_ctrl_idx]
        with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            this_res = pool.starmap(findCor, zip(range(0, len(gene_labels)), repeat(clust_bhve_idx, len(gene_labels))))
            clust15_res["BHVE_" + str(clust15)] = np.array(this_res)
        with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            this_res = pool.starmap(findCor, zip(range(0, len(gene_labels)), repeat(clust_ctrl_idx, len(gene_labels))))
            clust15_res["CTRL_" + str(clust15)] = np.array(this_res)
        clust15_dif_res[str(clust15)] = clust15_res["BHVE_" + str(clust15)] - clust15_res["CTRL_" + str(clust15)]

    clust15_df = pandas.DataFrame(clust15_dif_res, index=gene_labels)

    # 53 cluster level
    clust53_res = {}
    clust53_dif_res = {}
    # for clust53 in range(0, 53):
    #     print(clust53)
    #     clust_idx = np.where(cluster53_labels == clust53)[0]
    #     clust_bhve_idx = pandas.Series(bhve_idx).isin(clust_idx)
    #     clust_bhve_idx = bhve_idx[clust_bhve_idx]
    #     clust_ctrl_idx = pandas.Series(ctrl_idx).isin(clust_idx)
    #     clust_ctrl_idx = ctrl_idx[clust_ctrl_idx]
    #     with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
    #         this_res = pool.starmap(findCor, zip(range(0, len(gene_labels)), repeat(clust_bhve_idx, len(gene_labels))))
    #         clust53_res["BHVE_" + str(clust53)] = np.array(this_res)
    #     with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
    #         this_res = pool.starmap(findCor, zip(range(0, len(gene_labels)), repeat(clust_ctrl_idx, len(gene_labels))))
    #         clust53_res["CTRL_" + str(clust53)] = np.array(this_res)
    #     clust53_dif_res[str(clust53)] = clust53_res["BHVE_" + str(clust53)] - clust53_res["CTRL_" + str(clust53)]
    clust53_df = pandas.DataFrame(clust53_dif_res, index=gene_labels)

    return bulk_df, clust15_df, clust53_df

def main():
    start_time = time.perf_counter()
    gene_list, num_perm = parseArgs()

    # Read Data
    global data_mat
    global gene_labels
    global cond_labels
    global score_genes
    global score
    global score_genes_idx
    global cluster15_labels
    global cluster53_labels
    data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_data_mat.npz")
    gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_rownames.csv").iloc[:,
                  1].to_numpy()
    cond_labels = pandas.read_csv(
        "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_real_cond_labels.csv").iloc[:, 1].to_numpy()
    cluster15_labels = pandas.read_csv(
        "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_seuratclusters15.csv").iloc[:, 0].to_numpy()
    cluster53_labels = pandas.read_csv(
        "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_seuratclusters53.csv").iloc[:, 0].to_numpy()

    # Read in list of genes to create score with
    if gene_list == "ieg":
        score_genes = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/ieg_like_fos_egr1_npas4_detected_011521.csv").iloc[:, 0].to_numpy()
    elif gene_list == "neurogen":
        score_genes = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/neurogen_genes_final_050621.csv").iloc[:, 0].to_numpy()
    elif gene_list == "prog":
        score_genes = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/progenitor_sox2_nes_coexpress_051721.csv").iloc[:, 0].to_numpy()
    else:
        print("Invalid input gene list. The options are ieg, neurogen, and prog. Terminating program.")
        return(False)

    # Find the index of genes used to create the score
    score_genes_idx = pandas.Series(gene_labels).isin(score_genes)
    score_genes_idx = score_genes_idx[score_genes_idx].index

    # Read in Real Data
    real_bulk_df    = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/bulk_real_" + gene_list + "_score_cor_bvc.csv", index_col = 0)
    real_clust15_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/clust15_real_" + gene_list + "_score_cor_bvc.csv", index_col = 0)
    real_clust53_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/clust53_real_" + gene_list + "_score_cor_bvc.csv", index_col = 0)

    # # Real bhve and ctrl indexes
    # real_bhve_idx = np.where(cond_labels == "BHVE")[0]
    # real_ctrl_idx = np.where(cond_labels == "CTRL")[0]

    # Store values of indexes of genes that a positive and negative difference
    bulk_greater_idx = real_bulk_df['Dif'] > 0
    bulk_smaller_idx = real_bulk_df['Dif'] < 0
    clust15_bool_idx = {}  # key is cluster, value is a list of 2 (1. greater_idx and 2. smaller_idx)
    clust53_bool_idx = {}  # key is cluster, value is a list of 2 (1. greater_idx and 2. smaller_idx)
    for i in range(0, 15):
        greater_idx = real_clust15_df[str(i)] > 0
        smaller_idx = real_clust15_df[str(i)] < 0
        clust15_bool_idx[i] = [greater_idx, smaller_idx]
    for i in range(0, 53):
        greater_idx = real_clust53_df[str(i)] > 0
        smaller_idx = real_clust53_df[str(i)] < 0
        clust53_bool_idx[i] = [greater_idx, smaller_idx]

    # Prepare DataFrames that store the number of permutations more extreme than the real
    perm_greater_bulk = pandas.DataFrame(0, index = real_bulk_df.index, columns = ['nMoreExtreme'])
    perm_greater_clust15 = pandas.DataFrame(0, index=real_clust15_df.index, columns = real_clust15_df.columns)
    perm_greater_clust53 = pandas.DataFrame(0, index=real_clust53_df.index, columns = real_clust53_df.columns)

    time_before_perm = time.perf_counter()
    print(f"Time before starting permutations: {time.perf_counter() - start_time:0.4f} seconds")

    # Do permutations
    for i in range(0, num_perm):
        this_perm_start_time = time.perf_counter()
        print("Permutation " + str(i))
        perm_label = myShuffle(cond_labels)
        perm_bhve_idx = np.where(perm_label == "BHVE")[0]
        perm_ctrl_idx = np.where(perm_label == "CTRL")[0]
        global bhve_idx
        global ctrl_idx
        bhve_idx = perm_bhve_idx
        ctrl_idx = perm_ctrl_idx
        perm_bulk, perm_clust15, perm_clust53 = singleRun()

        # Compare the permuted results to the real results to see if they are more extreme
        perm_greater_bulk.loc[bulk_greater_idx, 'nMoreExtreme'] += perm_bulk.loc[bulk_greater_idx]['Dif'] > real_bulk_df.loc[bulk_greater_idx]['Dif']
        perm_greater_bulk.loc[bulk_smaller_idx, 'nMoreExtreme'] += perm_bulk.loc[bulk_smaller_idx]['Dif'] < real_bulk_df.loc[bulk_smaller_idx]['Dif']
        for i in range(0, 53):
            clust53_bool_idx = clust53_bool_idx[i]
            clust53_greater_idx = clust53_bool_idx[0]
            clust53_smaller_idx = clust53_bool_idx[1]
            # perm_greater_clust53.loc[clust53_greater_idx, str(i)] += perm_clust53.loc[clust53_greater_idx, str(i)] > real_clust53_df.loc[clust53_greater_idx, str(i)]
            # perm_greater_clust53.loc[clust53_smaller_idx, str(i)] += perm_clust53.loc[clust53_smaller_idx, str(i)] < real_clust53_df.loc[clust53_smaller_idx, str(i)]
            if i <= 14:
                clust15_bool_idx = clust15_bool_idx[i]
                clust15_greater_idx = clust15_bool_idx[0]
                clust15_smaller_idx = clust15_bool_idx[1]
                # print(clust15_bool_idx[0:5])
                # print(clust15_greater_idx[(len(clust15_greater_idx) - 5):len(clust15_greater_idx)])
                print(perm_greater_clust15)
                print(perm_clust15)
                print(real_clust15_df)
                print(clust15_greater_idx)
                print(perm_greater_clust15.columns)
                print(perm_clust15.columns)
                print(real_clust15_df.columns)
                perm_greater_clust15.loc[clust15_greater_idx, str(i)] += perm_clust15.loc[clust15_greater_idx, str(i)] > real_clust15_df.loc[clust15_greater_idx, str(i)]
                perm_greater_clust15.loc[clust15_smaller_idx, str(i)] += perm_clust15.loc[clust15_smaller_idx, str(i)] < real_clust15_df.loc[clust15_smaller_idx, str(i)]

        print(f"Time to complete permutation {i:0.1f}: {time.perf_counter() - this_perm_start_time:0.4f} seconds")

    perm_greater_bulk.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/bulk_perm_" + str(num_perm) + "_" + gene_list + "_score_cor_bvc.csv")
    perm_greater_clust15.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/clust15_perm_" + str(num_perm) + "_" + gene_list + "_score_cor_bvc.csv")
    perm_greater_clust53.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/clust53_perm_" + str(num_perm) + "_" + gene_list + "_score_cor_bvc.csv")
    print("Done")
    # bulk_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/bulk_real_neurogen_score_cor_bvc.csv")
    # clust15_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/clust15_real_neurogen_score_cor_bvc.csv")
    # clust53_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/clust53_real_neurogen_score_cor_bvc.csv")


if __name__ == '__main__':
    main()
