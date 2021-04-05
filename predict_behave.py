import sklearn
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import LogisticRegression
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
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
import pickle
import scipy
import pandas
import numpy

# infile = open("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/counts_seurat.pickle",'rb')
# df = pickle.load(infile)
# infile.close()

pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/t_sample_avg.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/c_hit_data.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/bulk_data.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/bulk1235_order_data.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/c_hit_order_data.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust_avg.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust53_avg_counts.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/c_hit_cluster_order_data_small_order.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/c_hit_cluster_order_data_small.txt",sep="\s")
df = pd_df.to_numpy()
# df = df[:,range(0,8)]
df.shape

pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust53_avg_counts.txt",sep="\s")
c_53_features = pd_df.columns
c_53 = pd_df.to_numpy()
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust53_avg.txt",sep="\s")
d_53_features = pd_df.columns
d_53 = pd_df.to_numpy()
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust_avg_counts.txt",sep="\s")
c_15_features = pd_df.columns
c_15 = pd_df.to_numpy()
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust_avg.txt",sep="\s")
d_15_features = pd_df.columns
d_15 = pd_df.to_numpy()
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/t_sample_avg_counts.txt",sep="\s")
c_b_features = pd_df.columns
c_b = pd_df.to_numpy()
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/t_sample_avg.txt",sep="\s")
d_b_features = pd_df.columns
d_b = pd_df.to_numpy()


# Meta Features
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_subsample_cluster_pct.txt",sep="\s")
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_sample_cluster_pct.txt",sep="\s")
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_replicate_cluster_pct_cells.csv")
pd_df = pd_df.iloc[:,1:]  # if csv has rownames
d_cells_features = pd_df.columns
d_cells = pd_df.to_numpy()

with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/bulk_biased_genes.txt','r') as f:
    d_b_biased = numpy.array(f.read().splitlines())

with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/15_biased_genes.txt','r') as f:
    d_15_biased = numpy.array(f.read().splitlines())

with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/53_biased_genes.txt','r') as f:
    d_53_biased = numpy.array(f.read().splitlines())

d_15_biased = [ '"' + x + '"' for x in d_15_biased ]
d_53_biased = [ '"' + x + '"' for x in d_53_biased ]

# sel = VarianceThreshold(threshold=0.005)
# df = sel.fit_transform(df)
# df.shape
# df = SelectKBest(chi2, k=500).fit_transform(df, [1,1,1,1,1,-1,-1,-1,-1,-1])
# df.shape


# df = d_cells
# features = d_cells_features
# df = pd_df[["rsrp1", "LOC101485431", "APK84-gp10", "smarce1"]].to_numpy()
# df = pd_df[['rsrp1', 'APK84-gp10', 'smarce1', 'LOC101475628', 'lrch3', 'LOC101484687', 'LOC101475431', 'LOC101473766', 'LOC106674590', 'map1b', 'LOC101479871', 'fam69a', 'slc2a4rg', 'LOC101466332', 'LOC112435883', 'LOC101469954', 'LOC101480168', 'LOC101482251', 'LOC101466366', 'pomgnt2']].to_numpy()
# df = pd_df[['rsrp1', 'APK84-gp10', 'smarce1', 'LOC101475628', 'lrch3', 'LOC101484687', 'LOC101475431', 'LOC101473766', 'LOC106674590', 'map1b', 'LOC101479871', 'fam69a', 'slc2a4rg', 'LOC101466332', 'LOC112435883']].to_numpy()
# sample = ["b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5"]
# pair = { 0:5, 1:6, 2:7, 3:8, 4:9 }
pairs = ["1", "2", "3", "4", "5"]
coef_df = pandas.DataFrame()
good_i = []
good_i_avg_prob = []
results = {}
success = False
for i in range(1, 2):
    # if i % 5 == 0:
    #     print("i: " + str(i))
    for j in range(1, 52):
        print("***")
        print("j: " + str(j))
        print("***")
        # if j % 5 == 0:
        #     print("j: " + str(j))
        results[str(i) + "," + str(j)] = []
        all_test_prob = []
        all_test_score = []
        for pair in pairs:
            print("Predicted Pair: b" + pair + " and c" + pair)
            sub_pair = [ x for x in list(pd_df.index) if x[1:2] == pair ]  # the subsamples that are in the test pair
            not_sub_pair = [ x for x in list(pd_df.index) if x[1:2] != pair ]  # the subsamples that are not in the test pair
            # sub_pair = [ x for x in list(pd_df.index) if x[2:3] == pair ]  # the subsamples that are in the test pair
            # not_sub_pair = [ x for x in list(pd_df.index) if x[2:3] != pair ]  # the subsamples that are not in the test pair
            xtrain = pd_df.loc[not_sub_pair]
            xtest = pd_df.loc[sub_pair]
            # ytrain = numpy.multiply([x[1:2] == "b" for x in not_sub_pair], 1)  # behave is 1 and control is 0
            # ytest = numpy.multiply([x[1:2] == "b" for x in sub_pair], 1)  # behave is 1 and control is 0
            ytrain = numpy.multiply([x[0:1] == "b" for x in not_sub_pair], 1)  # behave is 1 and control is 0
            ytest = numpy.multiply([x[0:1] == "b" for x in sub_pair], 1)  # behave is 1 and control is 0
            rc = LogisticRegression(C=1)
            # Biased
            # biased_idx = numpy.argwhere(numpy.isin(features, d_b_biased[j])).flatten().tolist()
            # biased_idx = numpy.argwhere(numpy.isin(features, d_b_biased[[i, j]])).flatten().tolist()
            # xtrain = xtrain[0:8, biased_idx]
            # xtest = xtest[0:2, biased_idx]
            # xtrain = sklearn.preprocessing.normalize(xtrain)
            # xtest = sklearn.preprocessing.normalize(xtest)
            # rc = RidgeClassifier(alpha = 0.9, normalize = True)
            # rc = RandomForestClassifier(n_estimators = 50, max_depth=2, random_state=0)
            # My All bhve > ctrl or all bhve < ctrl Feature Selection
            # b_mins = numpy.ndarray.min(xtrain[0:4,], axis = 0)
            # b_maxs = numpy.ndarray.max(xtrain[0:4,], axis = 0)
            # c_mins = numpy.ndarray.min(xtrain[4:8,], axis = 0)
            # c_maxs = numpy.ndarray.max(xtrain[4:8,], axis = 0)
            # b_big_idx = numpy.argwhere(b_mins > c_maxs)
            # c_big_idx = numpy.argwhere(c_mins > b_maxs)
            # all_big_idx = b_big_idx.flatten().tolist() + c_big_idx.flatten().tolist()
            # xtrain = xtrain[0:8, all_big_idx]
            # xtest = xtest[0:2, all_big_idx]
            # this_features = features[all_big_idx]
            # print(len(this_features))
            # My Feature Selection
            # greater_df = numpy.array([xtrain[x, :] > xtrain[x+4,:] for x in range(0,4)])
            # greater_idx = numpy.argwhere(numpy.logical_and(greater_df[0, :], numpy.logical_and(greater_df[1, :], numpy.logical_and(greater_df[2, :], greater_df[3, :]))))
            # smaller_df = numpy.array([xtrain[x, :] < xtrain[x+4,:] for x in range(0,4)])
            # smaller_idx = numpy.argwhere(numpy.logical_and(smaller_df[0, :], numpy.logical_and(smaller_df[1, :], numpy.logical_and(smaller_df[2, :], smaller_df[3, :]))))
            # all_idx = greater_idx.flatten().tolist() + smaller_idx.flatten().tolist()
            # greater_df = numpy.array([xtrain.iloc[x, :] > xtrain.iloc[x + int(xtrain.shape[0] / 2), :] for x in range(0,int(xtrain.shape[0] / 2))])
            # greater_idx = numpy.argwhere(numpy.sum(greater_df, axis = 0) > 0.33 * xtrain.shape[0])
            # smaller_df = numpy.array([xtrain.iloc[x, :] < xtrain.iloc[x + int(xtrain.shape[0] / 2), :] for x in range(0, int(xtrain.shape[0] / 2))])
            # smaller_idx = numpy.argwhere(numpy.sum(smaller_df, axis=0) > 0.33 * xtrain.shape[0])
            # all_idx = greater_idx.flatten().tolist() + smaller_idx.flatten().tolist()
            # xtrain = xtrain.iloc[:, all_idx]
            # xtest = xtest.iloc[:, all_idx]
            # this_features = features[all_idx]
            # print(len(all_idx))
            # 3. Pick features where behave and control are the most non-overlapping
            # b_min = xtrain.iloc[0:int(xtrain.shape[0] / 2), :].min(axis=0)
            # b_max = xtrain.iloc[0:int(xtrain.shape[0] / 2), :].max(axis=0)
            # c_min = xtrain.iloc[int(xtrain.shape[0] / 2):xtrain.shape[0], :].min(axis=0)
            # c_max = xtrain.iloc[int(xtrain.shape[0] / 2):xtrain.shape[0], :].max(axis=0)
            # b_up = b_min - c_max
            # c_up = c_min - b_max
            # b_up = b_up[b_up != 0]
            # c_up = c_up[c_up != 0]
            # all_up_features = list(b_up.index[(-b_up).argsort()[0:i]])  # features that are up in bhve
            # all_up_features.extend(list(c_up.index[(-c_up).argsort()[0:i]]))  # features that are up in ctrl
            # xtrain = xtrain[all_up_features]
            # xtest = xtest[all_up_features]
            # 4. Do a t-test and select the most significant ones
            # res = [scipy.stats.ttest_ind(xtrain.iloc[0:int(xtrain.shape[0] / 2), x],
            #                              xtrain.iloc[int(xtrain.shape[0] / 2):xtrain.shape[0], x]).pvalue for x in
            #        range(0, xtrain.shape[1])]
            # res = numpy.array(res)
            # res_idx = (res).argsort()
            # xtrain = xtrain.iloc[:, res_idx[0:i]]
            # xtest = xtest.iloc[:, res_idx[0:i]]
            # print(xtrain.columns)
            # minmax = MinMaxScaler().fit(xtrain)
            # xtrain = minmax.transform(xtrain)
            # xtest = minmax.transform(xtest)
            # b_sum = xtrain.iloc[0:int(xtrain.shape[0] / 2), :].sum(axis=0)
            # c_sum = xtrain.iloc[int(xtrain.shape[0] / 2):xtrain.shape[0], :].sum(axis=0)
            # dif = abs(b_sum - c_sum)
            # dif_idx = (-dif).argsort()
            # xtrain = xtrain.iloc[:, dif_idx[0:i]]
            # xtest = xtest.iloc[:, dif_idx[0:i]]
            # this_features = this_features[dif_idx[range(0, i)]]
            # scaler = StandardScaler()
            # a = scaler.fit(xtrain)
            # xtrain = pandas.DataFrame(data=scaler.transform(xtrain), index = xtrain.index, columns = xtrain.columns)
            # xtest = pandas.DataFrame(data=scaler.transform(xtest), index = xtest.index, columns = xtest.columns)
            # a = scaler.fit(xtest)  # according to the internet I shouldn't do this step
            a = rc.fit(xtrain, ytrain)
            # importance = rc.feature_importances_
            importance = numpy.abs(rc.coef_)
            sort_idx = (-importance).argsort()  # argsort give the indexes which would sort the array
            xtrain = xtrain.iloc[:, sort_idx[0, 0:j]]
            xtest = xtest.iloc[:, sort_idx[0, 0:j]]
            # print(importance[0,sort_idx[0, range(0,5)]])
            # xtrain = xtrain[0:8, sort_idx[0:5]]
            # xtest = xtest[0:2, sort_idx[0:5]]
            # this_features = this_features[sort_idx[0, range(0,j)]]
            a = rc.fit(xtrain, ytrain)
            # rc.predict(xtest)
            rc.predict_proba(xtest)
            # sfm = SelectFromModel(rc, threshold=1e-3)
            # sfm.fit(xtrain, ytrain)
            # xtrain = sfm.transform(xtrain)
            # xtest = sfm.transform(xtest)
            # conf_score = rc.decision_function(xtrain)
            # print("Confidence Score on train set: ", str(conf_score))
            # conf_score = rc.decision_function(xtest)
            # print("Confidence Score on test set: ", str(conf_score))
            # train_score = rc.score(xtrain, ytrain)
            # print("Accuracy on train set: ", str(train_score))
            test_score = rc.score(xtest, ytest)
            print("Accuracy on test set (Pair " + pair + "): ", str(test_score))
            print(xtrain.shape)
            # if test_score < 0.66:
            # #     print("Total Performance: " + str(len(all_test_score)) + "/" + str(len(all_test_score)+1))
            #     results[str(i) + "," + str(j)] = [-1, -1, -1, -1, -1]
            #     break
            # else:
            if test_score >= 1:
                all_test_score.append(test_score)
            results[str(i) + "," + str(j)].append(test_score)
            # # this_prob = rc.predict_proba(xtest)
            # # all_test_prob.append(this_prob[0][1])  # probability that the behave sample is behave
            # # all_test_prob.append(this_prob[1][0])  # probability that the control sample is control
            # # coef_df[pair] = xtest.columns
            if len(all_test_score) == 5:
                good_i.append([i,j])
                # good_i.append(j)
            #     good_i_avg_prob.append(numpy.mean(all_test_prob))
            #     if not success:
                print("Found a winner!" + str(i) + ", " + str(j))
                    # success = True
            #         break


results2 = pandas.DataFrame(data = results)
# results2 = results2.append(numpy.sum(results2, axis = 0)/results2.shape[0], ignore_index = True)
overall_pct = [ (results2.iloc[0,x] * 8 + results2.iloc[1,x] * 8 + results2.iloc[2,x] * 8 + results2.iloc[3,x] * 6 + results2.iloc[4,x] * 8)/38 for x in range(0, results2.shape[1])]
results2.loc[len(results2)] = overall_pct
best_idx = (-results2.iloc[results2.shape[0]-1, :]).argsort()
results2.iloc[:, best_idx[0:4]]

# Subsamples
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_subsample_data.txt",sep="\s")
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_subsample_counts.txt",sep="\s")
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_subsample_counts_15.txt",sep="\s")
pd_df_orig = pd_df

# TODO: Change my filter back
minmax = MinMaxScaler().fit(pd_df)
pd_df_scale = pandas.DataFrame(data = minmax.transform(pd_df), index = pd_df.index, columns = pd_df.columns)
pd_df2 = (pd_df == 0) * 1
mostly_nonzero_entries = pd_df2.sum(axis=0)[pd_df2.sum(axis=0) < 2]  # 2 for bulk and 5 for cluster
pd_df_scale = pd_df_scale[list(mostly_nonzero_entries.index)]
pd_df = pd_df_scale

# Select certain clusters
b_clust = ['.4', '.6', '.7', '.8', '.9']
b_clust_ind = [x for x in pd_df_orig.columns if x[-3:-1] in b_clust]
pd_df = pd_df_orig[b_clust_ind]

# Parallelized ML
import multiprocessing
from itertools import repeat
results = {}
for i in range(1, 200):
    if i % 5 == 0:
        print("i: " + str(i))
    j_min = 1
    j_max = 40
    pool = multiprocessing.Pool(20)
    this_result = pool.starmap(singleRun, zip( repeat(i), range(j_min,j_max) ))
    pool.close()
    for j in range(j_min, j_max):
        results[str(i) + ", " + str(j)] = this_result[j-j_min]

# View the Results
results2 = pandas.DataFrame(data = results)
overall_pct = [ (results2.iloc[0,x] * 8 + results2.iloc[1,x] * 8 + results2.iloc[2,x] * 8 + results2.iloc[3,x] * 6 + results2.iloc[4,x] * 8)/38 for x in range(0, results2.shape[1])]
results2.loc[len(results2)] = overall_pct
best_idx = (-results2.iloc[results2.shape[0]-1, :]).argsort()
results2.iloc[:, best_idx[0:4]]

# Plot the Results
from matplotlib import pyplot as plt
data = overall_pct
plt.hist(data, bins=30, alpha=0.5)
plt.title('Histogram of Accuracy in Iterations')
plt.xlabel('Accuracy')
plt.ylabel('Number of Iterations')
plt.savefig('iter_hist.png')

def singleRun(i, j):
    """

    :param i: Number of Features from my filter to pick
    :param j: Number of Top Features for the Model to pick
    :return: A list of the test scores for the 5 different models
    """
    results = []
    pairs = ["1", "2", "3", "4", "5"]
    for pair in pairs:
        # Split Data into Train and Test
        sub_pair = [x for x in list(pd_df.index) if x[2:3] == pair]  # the subsamples that are in the test pair
        not_sub_pair = [x for x in list(pd_df.index) if x[2:3] != pair]  # the subsamples that are not in the test pair
        xtrain = pd_df.loc[not_sub_pair]
        xtest = pd_df.loc[sub_pair]
        ytrain = numpy.multiply([x[1:2] == "b" for x in not_sub_pair], 1)  # behave is 1 and control is 0
        ytest = numpy.multiply([x[1:2] == "b" for x in sub_pair], 1)  # behave is 1 and control is 0
        # My Feature Selection
        # 1. Pick features that are in a consistent direction across subsamples
        # greater_df = numpy.array([xtrain.iloc[x, :] > xtrain.iloc[x + int(xtrain.shape[0] / 2), :] for x in range(0, int(xtrain.shape[0] / 2))])
        # greater_idx = numpy.argwhere(numpy.sum(greater_df, axis=0) > 0.33 * xtrain.shape[0])
        # smaller_df = numpy.array([xtrain.iloc[x, :] < xtrain.iloc[x + int(xtrain.shape[0] / 2), :] for x in range(0, int(xtrain.shape[0] / 2))])
        # smaller_idx = numpy.argwhere(numpy.sum(smaller_df, axis=0) > 0.33 * xtrain.shape[0])
        # all_idx = greater_idx.flatten().tolist() + smaller_idx.flatten().tolist()
        # xtrain = xtrain.iloc[:, all_idx]
        # xtest = xtest.iloc[:, all_idx]
        # 2. Pick features that have the biggest difference between behave and control
        b_sum = xtrain.iloc[0:int(xtrain.shape[0] / 2), :].sum(axis=0)
        c_sum = xtrain.iloc[int(xtrain.shape[0] / 2):xtrain.shape[0], :].sum(axis=0)
        dif = abs(b_sum - c_sum)
        dif_idx = (-dif).argsort()
        xtrain = xtrain.iloc[:, dif_idx[0:i]]
        xtest = xtest.iloc[:, dif_idx[0:i]]
        # Standardize the Data
        # scaler = StandardScaler()
        # a = scaler.fit(xtrain)
        # scaler = MinMaxScaler().fit(xtrain)
        # xtrain = pandas.DataFrame(data=scaler.transform(xtrain), index=xtrain.index, columns=xtrain.columns)
        # xtest = pandas.DataFrame(data=scaler.transform(xtest), index=xtest.index, columns=xtest.columns)
        # Fit the Model
        rc = LogisticRegression(C=1)
        a = rc.fit(xtrain, ytrain)
        # Select the Top Features
        importance = numpy.abs(rc.coef_)
        sort_idx = (-importance).argsort()  # argsort give the indexes which would sort the array
        xtrain = xtrain.iloc[:, sort_idx[0, 0:j]]
        xtest = xtest.iloc[:, sort_idx[0, 0:j]]
        # Refit the model
        a = rc.fit(xtrain, ytrain)
        test_score = rc.score(xtest, ytest)
        # if test_score < 0.66:
        #     results = [-1, -1, -1, -1, -1]
        #     break
        results.append(test_score)
    return results


d_15_biased[numpy.argmax(good_i_avg_prob)]
numpy.argwhere(numpy.isin(features, d_15_biased[j])).flatten().tolist()

        # print(rc.n_features_)
        # print(len(rc.coef_[0,:]))
        # coef_df[b_pair] = rc.coef_[0,:]
        # coef_df[b_pair] = pd_df.columns[sort_idx[0, range(0, i)]]
        # coef_df[b_pair] = this_features.tolist()




import csv
coef_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/coef_cluster_unbiased_44.txt", quoting=csv.QUOTE_NONE)
coef_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/coef_bulk_unbiased_299.txt", quoting=csv.QUOTE_NONE)
coef_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/coef_cluster53_unbiased_358.txt", quoting=csv.QUOTE_NONE)
coef_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/coef_cluster15_unbiased_logist.txt", quoting=csv.QUOTE_NONE)
coef_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/coef_cluster15_unbiased_logist_27.txt", quoting=csv.QUOTE_NONE)
coef_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/coef_bulk_3_172021.txt", quoting=csv.QUOTE_NONE)
coef_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/coef_15_2_172021.txt", quoting=csv.QUOTE_NONE)
coef_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/coef_53_4_172021.txt", quoting=csv.QUOTE_NONE)

# cv_scores = cross_val_score(rc, xtrain, ytrain, cv=4)
# print("CV average score: %.2f" % cv_scores.mean())
# cm = confusion_matrix(ytest, ypred)
# print(cm)

importance = numpy.abs(rc.coef_)
sort_idx = (-importance).argsort()
importance[0, sort_idx]

# xtrain = numpy.array([test[0], test[1], test[2], test[3], test[5], test[6], test[7], test[8]])
# xtest = numpy.array([test[4], test[9]])
# ytrain = numpy.array([1, 1, 1, 1, 0, 0, 0, 0])
# ytest = numpy.array([1, 0])

rc = RidgeClassifier()
rc.fit(xtrain, ytrain)
score = rc.score(xtrain, ytrain)
print("Score: ", score)

cv_scores = cross_val_score(rc, xtrain, ytrain, cv=4)
print("CV average score: %.2f" % cv_scores.mean())

ypred = rc.predict(xtest)

cm = confusion_matrix(ytest, ypred)
print(cm)

cr = classification_report(ytest, ypred)
print(cr)


# sample = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample.txt",sep="\s")
# sample = sample.transpose()

# b1 = df[sample == "b1"]
# b2 = df[sample == "b2"]
# b3 = df[sample == "b3"]
# b4 = df[sample == "b4"]
# b5 = df[sample == "b5"]
# c1 = df[sample == "c1"]
# c2 = df[sample == "c2"]
# c3 = df[sample == "c3"]
# c4 = df[sample == "c4"]
# c5 = df[sample == "c5"]
# test = numpy.array(b1, b2)

# outfile = open("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/np_counts.pickle",'wb')
# pickle.dump(test, outfile)
# outfile.close()

this_pair = [3, 8]
xtrain = [2,7]
xtrain = df[xtrain]
xtest = df[this_pair]
ytrain = numpy.array([1,-1])
ytest = numpy.array([1, -1])
rc = RidgeClassifier(alpha = 1, normalize = False)
rc.fit(xtrain, ytrain)
conf_score = rc.decision_function(xtest)
print("Confidence Score on test set: ", str(conf_score))
rc.coef_[0,1:10]


"""
For Cluster Independent Analysis
"""
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust53_avg.txt",sep="\s")
max_iter = 1000
for clust in range(15,53):
    print("Cluster " + str(clust))
    success = False
    df = pd_df.iloc[:, (0+32471*clust):(32471*(clust+1))].to_numpy()
    sample = ["b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5"]
    pair = { 0:5, 1:6, 2:7, 3:8, 4:9 }
    coef_df = pandas.DataFrame()
    for i in range(1, max_iter):
        all_train_score = []
        all_test_score = []
        for b_pair, c_pair in pair.items():
            # print("Predicted Pair: " + sample[b_pair] + " and " + sample[c_pair])
            this_pair = [b_pair, c_pair]
            # Define the training and test sets
            xtrain = [ x for x in range(0,10) if x not in this_pair ]
            xtrain = df[xtrain]
            xtest = df[this_pair]
            ytrain = numpy.array([1, 1, 1, 1, -1, -1, -1, -1])
            ytest = numpy.array([1, -1])
            # Standardize the features
            scaler = StandardScaler()
            a = scaler.fit(xtrain)
            xtrain = scaler.transform(xtrain)
            a = scaler.fit(xtest)
            xtest = scaler.transform(xtest)
            # Create a linear model and train it
            rc = LogisticRegression()
            a = rc.fit(xtrain, ytrain)
            # Find the i most important features
            importance = numpy.abs(rc.coef_)
            sort_idx = (-importance).argsort()
            xtrain = xtrain[0:8, sort_idx[0, range(0,i)]]
            xtest = xtest[0:8, sort_idx[0, range(0,i)]]
            # Retrain the model and find the test accuracy
            a = rc.fit(xtrain, ytrain)
            test_score = rc.score(xtest, ytest)
            if test_score == 1:
                all_test_score.append(test_score)
        if len(all_test_score) == 5:
            print("Minimum Number of Genes: " + str(i))
            success = True
            break
    if not success:
        print("None of the Top 0-" + str(max_iter) + " Genes Correctly Predict All Samples")

