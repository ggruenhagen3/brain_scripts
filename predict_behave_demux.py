# Load Libraries
import sklearn
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import BaggingRegressor
from sklearn import linear_model
from sklearn.datasets import load_iris
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_regression
from sklearn.feature_selection import mutual_info_regression
from sklearn.feature_selection import SelectKBest
from statsmodels.sandbox.stats.multicomp import multipletests
import pickle
import scipy
from itertools import repeat
from itertools import product
import multiprocessing
import pandas
import numpy
import sklearn
from matplotlib import pyplot as plt
import seaborn as sns
pandas.options.mode.chained_assignment = None  # default='warn'

# Load Data
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/demux_data/cluster0_data.csv", index_col = 0)
sub_meta = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/demux_data/subsample_meta2.csv", index_col = 0)
sub_meta['cond_bin'] = numpy.array(sub_meta['cond'] == 'BHVE') * 1

"""
Individual =====================================================================================================
"""
# Predict 1 Individual at a Time
params = {'bhve_metric': ['bower_activity_index'],
          'num_biggest': [0, 20, 40, 45, 50, 55, 60, 65, 70, 75],
          'reg_type':    ['lin', 'lasso', 'bag'],
          'scale_type': ['none', 'std', 'minmax'],
          'num_k': [0, 10, 20, 25, 50],
          'doPCA': [False]}
big_df = pandas.DataFrame(list(product(params['bhve_metric'], params['num_biggest'], params['reg_type'], params['scale_type'], params['num_k'], params['doPCA'])),
                          columns=['bhve_metric', 'num_biggest', 'reg_type', 'scale_type', 'num_k', 'doPCA'])
bad_idx = big_df.loc[(big_df['num_biggest'] != 0) & (big_df['num_k']  > big_df['num_biggest'])].index
big_df = big_df.iloc[~ big_df.index.isin(bad_idx) ].reset_index(drop=True)
big_df[["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]] = 0
for i in big_df.index:
    print(str(i) + "/" + str(len(big_df.index)))
    bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA = big_df.loc[i, ["bhve_metric", "num_biggest", "reg_type", 'scale_type', 'num_k', 'doPCA']]
    # bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA = ['cond_bin', 50, 'rc', 'none', 0, True]
    all_res = sub_meta[['subsample', 'cond']]
    all_res['real'] = sub_meta[bhve_metric]
    all_res['pred'] = 0
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        all_ypred = pool.starmap(singleIndPred, zip(pd_df.index, repeat(bhve_metric), repeat(num_biggest), repeat(reg_type), repeat(scale_type), repeat(num_k), repeat(doPCA)))
        all_res['pred'] = numpy.concatenate(all_ypred)
        big_df.loc[i, ["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]] = analyzeRes(all_res)


big_df.sort_values('all_cor', ascending = False)

"""
Pair =====================================================================================================
"""
# Predict 1 Pair at a Time
all_pairs = list(set(sub_meta['pair']))
all_pair_idx = numpy.concatenate([ sub_meta.loc[sub_meta['pair'].isin([p]),].index for p in all_pairs ])
params = {'bhve_metric': ['bower_activity_index'],
          'num_biggest': [30, 35, 45, 55, 60, 65, 70, 75],
          'reg_type':    ['lin', 'lasso', 'bag'],
          'scale_type': ['none', 'std', 'minmax'],
          'num_k': [0, 25, 30, 35, 40, 50],
          'doPCA': [False]}
big_df = pandas.DataFrame(list(product(params['bhve_metric'], params['num_biggest'], params['reg_type'], params['scale_type'], params['num_k'], params['doPCA'])),
                          columns=['bhve_metric', 'num_biggest', 'reg_type', 'scale_type', 'num_k', 'doPCA'])
bad_idx = big_df.loc[(big_df['num_biggest'] != 0) & (big_df['num_k']  > big_df['num_biggest'])].index
big_df = big_df.iloc[~ big_df.index.isin(bad_idx) ].reset_index(drop=True)
big_df.loc[big_df['reg_type'].isin(['logreg', 'rc', 'svm', 'rf', 'gaus', 'tree', 'knn', 'lasso_bin']), 'bhve_metric'] = 'cond_bin'
big_df[["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]] = 0
for i in big_df.index:
    print(str(i) + "/" + str(len(big_df.index)))
    bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA = big_df.loc[i, ["bhve_metric", "num_biggest", "reg_type", 'scale_type', 'num_k', 'doPCA']]
    # bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA = big_df.loc[217, ["bhve_metric", "num_biggest", "reg_type", 'scale_type', 'num_k', 'doPCA']]
    # bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA = ['bower_activity_index', 45, 'bag', 'std', 0, False]
    all_res = sub_meta[['subsample', 'cond']]
    all_res['real'] = sub_meta[bhve_metric]
    all_res['pred'] = 0
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        all_ypred = pool.starmap(pairPred, zip(all_pairs, repeat(bhve_metric), repeat(num_biggest), repeat(reg_type), repeat(scale_type), repeat(num_k), repeat(doPCA)))
        all_res.loc[all_pair_idx, 'pred'] = numpy.concatenate(all_ypred)
        big_df.loc[i, ["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]] = analyzeRes(all_res)


best_idx = (-big_df['all_cor']).argsort()
big_df.loc[best_idx[0:10], ["bhve_metric", "num_biggest", "reg_type", 'scale_type', 'num_k', 'doPCA']]
big_df.loc[best_idx[0:10], ["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]]


big_df.to_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/pred_big_df.csv')
test = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/pred_big_df.csv", index_col = 0)
test_idx = (-test['all_cor']).argsort()
test.loc[test_idx[0:4], ["bhve_metric", "num_biggest", "reg_type", 'scale_type', 'num_k', 'doPCA']]
test.loc[test_idx[0:4], ["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]]

plotRes(all_res)

"""
Pair Final =================================================================================================
"""
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/demux_data/cluster0_data.csv", index_col = 0)
all_pairs = list(set(sub_meta['pair']))
all_pair_idx = numpy.concatenate([ sub_meta.loc[sub_meta['pair'].isin([p]),].index for p in all_pairs ])
all_pairs_genes = pandas.DataFrame(data={'pair': all_pairs})
all_res = sub_meta[['subsample', 'cond']]
all_res['real'] = sub_meta['cond_bin']
all_res['pred'] = 0
all_res['conf'] = 0
all_ypred = []
all_conf = []
for x in range(1,51):
    all_pairs_genes['gene' + str(x)] = 0


for p in all_pairs:
    train_loc = sub_meta.loc[~sub_meta['pair'].isin([p]),].index
    test_loc = sub_meta.loc[sub_meta['pair'].isin([p]),].index
    xtrain = pd_df.loc[train_loc, ]
    xtest  = pd_df.loc[test_loc, ]
    ytrain = sub_meta.loc[train_loc, 'cond_bin']
    ytest = sub_meta.loc[test_loc, 'cond_bin']
    xtrain, xtest = biggestBVCDif(xtrain, xtest, num_biggest=50)
    rc = RidgeClassifier().fit(xtrain, ytrain)
    all_ypred.append(rc.predict(xtest))
    all_conf.append(rc.decision_function(xtest))
    this_genes = xtest.columns
    all_pairs_genes.loc[all_pairs_genes['pair'] == p, all_pairs_genes.columns[1:51]] = this_genes


all_res.loc[all_pair_idx, 'pred'] = numpy.concatenate(all_ypred)
all_res.loc[all_pair_idx, 'conf'] = numpy.concatenate(all_conf)
analyzeRes(all_res)
# all_pairs_genes.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/best_pred_clust0_genes.csv")

"""
Pair RC =====================================================================================================
"""
all_pairs = list(set(sub_meta['pair']))
all_pair_idx = numpy.concatenate([ sub_meta.loc[sub_meta['pair'].isin([p]),].index for p in all_pairs ])
params = {'num_biggest': [0],
          'num_biggest2': [40, 45, 50, 55, 60, 65, 70, 75],
          'scale_type': ['none', 'std', 'minmax'],
          'num_k2': [0, 5, 10, 15, 20, 25, 40, 50],
          'my_a': [0.5, 0.75, 1, 1.25, 1.5]}
big_df = pandas.DataFrame(list(product(params['num_biggest'], params['num_biggest2'], params['scale_type'], params['num_k2'], params['my_a'])),
                          columns=['num_biggest', 'num_biggest2', 'scale_type', 'num_k2', 'my_a'])
bad_idx = big_df.loc[(big_df['num_biggest'] != 0) & (big_df['num_k2'] > big_df['num_biggest'])].index
bad_idx = big_df.loc[(big_df['num_biggest2'] != 0) & (big_df['num_k2'] > big_df['num_biggest2'])].index
big_df = big_df.iloc[~ big_df.index.isin(bad_idx) ].reset_index(drop=True)
big_df[["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]] = 0
for i in big_df.index:
    print(str(i) + "/" + str(len(big_df.index)))
    num_biggest, num_biggest2, scale_type, num_k2, my_a = big_df.loc[i, ["num_biggest", 'num_biggest2', 'scale_type', 'num_k2', 'my_a']]
    # num_biggest, scale_type, num_k2, my_a = [50, 'none', 5, 1]
    all_res = sub_meta[['subsample', 'cond']]
    all_res['real'] = sub_meta[bhve_metric]
    all_res['pred'] = 0
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        all_ypred = pool.starmap(pairPredRC, zip(all_pairs, repeat(num_biggest), repeat(num_biggest2), repeat(scale_type), repeat(num_k2), repeat(my_a)))
        all_res.loc[all_pair_idx, 'pred'] = numpy.concatenate(all_ypred)
        big_df.loc[i, ["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]] = analyzeRes(all_res)


big_df.sort_values('all_cor', ascending = False)

"""
Pool =====================================================================================================
"""
# Predict 1 Pool at a Time
all_pools = list(set(sub_meta['pool']))
all_pools_idx = numpy.concatenate([ sub_meta.loc[sub_meta['pool'].isin([p]),].index for p in all_pools ])
params = {'bhve_metric': ['cond_bin'],
          'num_biggest': [0, 20, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100],
          'reg_type':    ['rc', 'logreg', 'logregl1', 'rf', 'tree', 'knn'],
          'scale_type': ['none', 'std', 'minmax'],
          'num_k': [0, 10, 20, 25, 50],
          'doPCA': [False]}
big_df = pandas.DataFrame(list(product(params['bhve_metric'], params['num_biggest'], params['reg_type'], params['scale_type'], params['num_k'], params['doPCA'])),
                          columns=['bhve_metric', 'num_biggest', 'reg_type', 'scale_type', 'num_k', 'doPCA'])
bad_idx = big_df.loc[(big_df['num_biggest'] != 0) & (big_df['num_k']  > big_df['num_biggest'])].index
big_df = big_df.iloc[~ big_df.index.isin(bad_idx) ].reset_index(drop=True)
big_df[["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]] = 0
for i in big_df.index:
    print(str(i) + "/" + str(len(big_df.index)))
    bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA = big_df.loc[i, ["bhve_metric", "num_biggest", "reg_type", 'scale_type', 'num_k', 'doPCA']]
    # bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA = ['depth', 0, 'lasso', 'none', 0, False]
    all_res = sub_meta[['subsample', 'cond']]
    all_res['real'] = sub_meta[bhve_metric]
    all_res['pred'] = 0
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        all_ypred = pool.starmap(poolPred, zip(all_pools, repeat(bhve_metric), repeat(num_biggest), repeat(reg_type), repeat(scale_type), repeat(num_k), repeat(doPCA)))
        all_res.loc[all_pools_idx, 'pred'] = numpy.concatenate(all_ypred)
        big_df.loc[i, ["mean_b", "mean_c", "mean_dif", "num_b_in_c", "num_c_in_b", "all_cor"]] = analyzeRes(all_res)

big_df.sort_values('all_cor', ascending = False)


# for bhve_metric in params['bhve_metric']:
#     for num_biggest in params['num_biggest']:
#         for reg_type in params['reg_type']:

# for this_sub in pd_df.index:
#     print(this_sub)
#     all_res.loc[this_sub, 'pred'] = singleIndPred(this_sub, "depth", 1000, "lasso")
#
#
# all_res['real_pred_dif'] = all_res['pred'] - all_res['depth']
# printResStats(all_res)

def plotRes(all_res, fname1 = "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/pred_box_1.png", fname2 = "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/pred_scatter_1.png"):
    # Box Plot
    sns.set(style="whitegrid")
    plt.figure()
    sns.boxplot(x="cond", y="pred", data=all_res, showfliers=False)
    sns.swarmplot(x="cond", y="pred", data=all_res, color=".25")
    plt.show()
    plt.savefig(fname1)
    plt.cla()
    plt.clf()
    # Scatter Plot
    sns.set(style="whitegrid")
    plt.figure()
    sns.scatterplot(data=all_res, x="real", y="pred", hue="cond")
    plt.show()
    plt.savefig(fname2)
    plt.cla()
    plt.clf()


def poolPred(pool_pair, bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA):
    train_loc = sub_meta.loc[~sub_meta['pool'].isin([pool_pair]),].index
    test_loc = sub_meta.loc[sub_meta['pool'].isin([pool_pair]),].index
    xtrain = pd_df.loc[train_loc, ]
    xtest  = pd_df.loc[test_loc, ]
    ytrain = sub_meta.loc[train_loc, bhve_metric]
    ytest = sub_meta.loc[test_loc, bhve_metric]
    if num_biggest != 0:
        xtrain, xtest = biggestBVCDif(xtrain, xtest, num_biggest=num_biggest)
    if doPCA:
        xtrain, xtest = myPCA(xtrain, xtest)
    if num_k != 0:
        if num_k == 'p':
            xtrain, xtest = myFeatSel(xtrain, xtest, ytrain, type=num_k)
        else:
            xtrain, xtest = myFeatSel(xtrain, xtest, ytrain, num_k=num_k)
    if scale_type != "none":
        xtrain, xtest = myScale(xtrain, xtest, type=scale_type)
    ypred = myReg(xtrain, xtest, ytrain, ytest, type=reg_type)
    return ypred


def pairPred(pair_num, bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA):
    train_loc = sub_meta.loc[~sub_meta['pair'].isin([pair_num]),].index
    test_loc = sub_meta.loc[sub_meta['pair'].isin([pair_num]),].index
    xtrain = pd_df.loc[train_loc, ]
    xtest  = pd_df.loc[test_loc, ]
    ytrain = sub_meta.loc[train_loc, bhve_metric]
    ytest = sub_meta.loc[test_loc, bhve_metric]
    if num_biggest != 0:
        xtrain, xtest = biggestBVCDif(xtrain, xtest, num_biggest=num_biggest)
    if doPCA:
        xtrain, xtest = myPCA(xtrain, xtest)
    if num_k != 0:
        if num_k == 'p':
            xtrain, xtest = myFeatSel(xtrain, xtest, ytrain, type=num_k)
        else:
            xtrain, xtest = myFeatSel(xtrain, xtest, ytrain, num_k=num_k)
    if scale_type != "none":
        xtrain, xtest = myScale(xtrain, xtest, type=scale_type)
    ypred = myReg(xtrain, xtest, ytrain, ytest, type=reg_type)
    return ypred


def pairPredRC(pair_num, num_biggest, num_biggest2, scale_type, num_k2, my_a = 1, my_s = 'auto'):
    train_loc = sub_meta.loc[~sub_meta['pair'].isin([pair_num]),].index
    test_loc = sub_meta.loc[sub_meta['pair'].isin([pair_num]),].index
    xtrain = pd_df.loc[train_loc, ]
    xtest  = pd_df.loc[test_loc, ]
    ytrain = sub_meta.loc[train_loc, bhve_metric]
    ytest = sub_meta.loc[test_loc, bhve_metric]
    if num_biggest != 0:
        xtrain, xtest = biggestBVCDif(xtrain, xtest, num_biggest=num_biggest)
    if num_biggest2 != 0:
        xtrain, xtest = biggestBVCDif2(xtrain, xtest, num_biggest2=num_biggest2)
    if scale_type != "none":
        xtrain, xtest = myScale(xtrain, xtest, type=scale_type)
    rc = RidgeClassifier(alpha=my_a, solver=my_s).fit(xtrain, ytrain)
    if num_k2 != 0:
        imp = abs(rc.coef_)
        xtrain = xtrain.iloc[0::, (-imp).argsort()[0][0:num_k2]]
        xtest = xtest.iloc[0::, (-imp).argsort()[0][0:num_k2]]
    rc = rc.fit(xtrain, ytrain)
    ypred = rc.predict(xtest)
    # print(rc.predict(xtest))
    # print(rc.decision_function(xtest))
    return ypred


def singleIndPred(this_sub, bhve_metric, num_biggest, reg_type, scale_type, num_k, doPCA):
    train_loc = pd_df.index[pd_df.index != this_sub]
    test_loc  = pd_df.index[pd_df.index == this_sub]
    xtrain = pd_df.loc[train_loc,]
    xtest  = pd_df.loc[test_loc,]
    ytrain = sub_meta.loc[train_loc, bhve_metric]
    ytest  = sub_meta.loc[test_loc, bhve_metric]
    if num_biggest != 0:
        xtrain, xtest = biggestBVCDif(xtrain, xtest, num_biggest = num_biggest)
    if doPCA:
        xtrain, xtest = myPCA(xtrain, xtest)
    if num_k != 0:
        if num_k == 'p':
            xtrain, xtest = myFeatSel(xtrain, xtest, ytrain, type = num_k)
        else:
            xtrain, xtest = myFeatSel(xtrain, xtest, ytrain, num_k = num_k)
    if scale_type != "none":
        xtrain, xtest = myScale(xtrain, xtest, type = scale_type)
    ypred = myReg(xtrain, xtest, ytrain, ytest, type = reg_type)
    return ypred


def analyzeRes(all_res):
    b_names = all_res.index[numpy.char.startswith(all_res.index.tolist(), 'b')]
    c_names = all_res.index[numpy.char.startswith(all_res.index.tolist(), 'c')]
    if len(set(all_res['pred'])) == 2:
        # How to Analyze Results for Binary Classification
        mean_b = 0
        mean_c = 0
        mean_dif = 0
        # num_b_in_c = len(all_res.loc[(all_res['real'] == 'CTRL') & (all_res['pred'] == 'BHVE'),].index)
        # num_c_in_b = len(all_res.loc[(all_res['real'] == 'BHVE') & (all_res['pred'] == 'CTRL'),].index)
        num_b_in_c = len(all_res.loc[(all_res['real'] == 0) & (all_res['pred'] == 1),].index)
        num_c_in_b = len(all_res.loc[(all_res['real'] == 1) & (all_res['pred'] == 0),].index)
        all_cor = 1 - (num_b_in_c + num_c_in_b)/len(all_res.index)
    else:
        # How to Analyze Results for Continuous Classification
        mean_b = numpy.mean(all_res.loc[b_names, 'pred'])
        mean_c = numpy.mean(all_res.loc[c_names, 'pred'])
        mean_dif = mean_b - mean_c
        b_min = numpy.min(all_res.loc[b_names, 'pred'])
        b_max = numpy.max(all_res.loc[b_names, 'pred'])
        c_min = numpy.min(all_res.loc[c_names, 'pred'])
        c_max = numpy.max(all_res.loc[c_names, 'pred'])
        num_b_in_c = len(all_res.loc[b_names, 'pred'][all_res.loc[b_names, 'pred'] < c_max])
        num_c_in_b = len(all_res.loc[c_names, 'pred'][all_res.loc[c_names, 'pred'] > b_min])
        all_cor = numpy.corrcoef(all_res['pred'], all_res['real'])[0, 1]
        b_cor = numpy.corrcoef(all_res.loc[b_names, 'pred'], all_res.loc[b_names, 'real'])[0, 1]
    return mean_b, mean_c, mean_dif, num_b_in_c, num_c_in_b, all_cor


def printResStats(all_res):
    all_b_sort = numpy.sort(all_res.loc[all_res.index[numpy.char.startswith(all_res.index.tolist(), 'b')], 'pred'])
    all_c_sort = numpy.sort(all_res.loc[all_res.index[numpy.char.startswith(all_res.index.tolist(), 'c')], 'pred'])
    print("B Sort:")
    print(all_b_sort)
    print("C Sort:")
    print(all_c_sort)
    print("BHVE Mean: " + str(numpy.mean(all_b_sort)) + ", CTRL Mean: " + str(numpy.mean(all_c_sort)))


def myScale(xtrain, xtest, type = "std"):
    if type == "std":
        scaler = StandardScaler().fit(xtrain)
    elif type == "minmax":
        scaler = MinMaxScaler().fit(xtrain)
    else:
        print("Not an valid scaler type. Please select std or minmax.")
    xtrain = pandas.DataFrame(data=scaler.transform(xtrain), index=xtrain.index, columns=xtrain.columns)
    xtest = pandas.DataFrame(data=scaler.transform(xtest), index=xtest.index, columns=xtest.columns)
    return xtrain, xtest

def myFeatSel(xtrain, xtest, ytrain, num_k = 50, type = "k", model = "f_reg"):
    # Model for Feature Selection
    this_model = f_regression
    if model == "f_reg":
        this_model = f_regression
    elif model == "mutual":
        this_model = mutual_info_regression
    else:
        print("Not a valid model. Please select f_reg or mutual.")
    # Feature Selection type
    sel = SelectKBest(this_model, k=2)
    sel = sel.fit(xtrain, ytrain)
    if type == "k":
        gene_scores_sort = (-sel.scores_).argsort()
        xtrain = xtrain.iloc[:, gene_scores_sort[0:num_k]]
        xtest  = xtest.iloc[:, gene_scores_sort[0:num_k]]
    elif type == "p":
        xtrain = xtrain.iloc[:, numpy.where(sel.pvalues_ < 0.05)[0]]
        xtest  = xtest.iloc[:, numpy.where(sel.pvalues_ < 0.05)[0]]
    elif type == "fdr":  # no significant fdr genes :(
        qvalues = multipletests(pvals=sel.pvalues_, method="fdr_bh")
        qvalues = multipletests(pvals=sel.pvalues_[numpy.logical_not(numpy.isnan(sel.pvalues_))], method="fdr_bh")
    else:
        print("Not a valid type. Please select k or p. FDR is possible, but there's no sig genes, so I haven't implemented it.")
    if ((type == "p" or type == "fdr") and model == "mutual"):
        print("fdr and mutual cannot be used at once because mutual does not produce pvalues.")
    return xtrain, xtest

def myPCA(xtrain, xtest):
    pca_m = PCA()
    pca_m = pca_m.fit(xtrain)
    xtrain = pandas.DataFrame(data=pca_m.transform(xtrain), index = xtrain.index)
    xtest = pandas.DataFrame(data=pca_m.transform(xtest), index = xtest.index)
    return xtrain, xtest


def biggestBVCDif(xtrain, xtest, num_biggest = 50):
    b_sum = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'b')],].sum(axis=0)
    c_sum = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'c')],].sum(axis=0)
    dif = abs(b_sum - c_sum)
    dif_idx = (-dif).argsort()
    xtrain = xtrain.iloc[:, dif_idx[0:num_biggest]]
    xtest = xtest.iloc[:, dif_idx[0:num_biggest]]
    return xtrain, xtest

def biggestBVCDif2(xtrain, xtest, num_biggest2 = 50):
    b_median = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'b')],].median(axis=0)
    c_median = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'c')],].median(axis=0)
    # b_mean = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'b')],].mean(axis=0)
    # c_mean = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'c')],].mean(axis=0)
    # b_max = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'b')],].max(axis=0)
    # c_max = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'c')],].max(axis=0)
    # b_min = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'b')],].min(axis=0)
    # c_min = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'c')],].min(axis=0)
    # bmin_g_cmax = b_min > c_max
    # # Case 1
    # b_g_c_mean = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'b')],] > pandas.Series(c_mean, index=xtrain.columns)
    # b_l_c_max = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'b')],] < pandas.Series(c_max, index=xtrain.columns)
    # b_bad_count = b_l_c_max & b_g_c_mean
    # b_bad_count = b_bad_count.sum(axis = 0)
    # # Case 2
    # c_g_b_mean = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'c')],] > pandas.Series(b_mean, index=xtrain.columns)
    # c_l_b_max = xtrain.loc[xtrain.index[numpy.char.startswith(xtrain.index.tolist(), 'c')],] < pandas.Series(b_max, index=xtrain.columns)
    # c_bad_count = c_l_b_max & c_g_b_mean
    # c_bad_count = c_bad_count.sum(axis = 0)
    # Try 3
    dif = abs(b_median - c_median)
    dif_idx = (-dif).argsort()
    xtrain = xtrain.iloc[:, dif_idx[0:num_biggest2]]
    xtest = xtest.iloc[:, dif_idx[0:num_biggest2]]
    return xtrain, xtest

def myReg(xtrain, xtest, ytrain, ytest, type = "lin"):
    # Fit Model
    if type == "lin":
        rc = LinearRegression()
    elif type == "dtr":
        rc = DecisionTreeRegressor()  # not continuous
    elif type == "bag":
        rc = BaggingRegressor()
    elif type == "lasso":
        rc = Lasso()
        # Last Continuous Classfier #
    elif type == "logreg":
        rc = LogisticRegression(C=1, solver = 'liblinear')
    elif type == "rc":
        rc = RidgeClassifier(alpha = 1, solver = 'auto')
    elif type == "svm":
        rc = svm.LinearSVC()
    elif type == "rf":
        rc = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
    elif type == "gaus":
        rc = GaussianNB()
    elif type == "knn":
        rc = KNeighborsClassifier()
    elif type == "tree":
        rc = DecisionTreeClassifier()
    elif type == "logregl1":
        # Aka LASSO
        rc = LogisticRegression(C=1, solver='liblinear', penalty='l1')
    else:
        print("Not an valid regression type. Valid Binary Classifiers are: logreg, rc, svm, rf, gaus, knn, logregl1, or tree. Valid Continuous Classifiers are: lin, bag, lasso or dtr.")
    a = rc.fit(xtrain, ytrain)
    ypred = rc.predict(xtest)
    # ypred = rc.predict(xtest)[0]
    return ypred

