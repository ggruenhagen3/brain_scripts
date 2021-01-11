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
import pandas
import numpy

# infile = open("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/counts_seurat.pickle",'rb')
# df = pickle.load(infile)
# infile.close()

# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/t_sample_avg.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/c_hit_data.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/bulk_data.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/bulk1235_order_data.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/c_hit_order_data.txt",sep="\s")
# pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust_avg.txt",sep="\s")
pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/grn/bb_info/sample_clust53_avg_counts.txt",sep="\s")
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

# df = pd_df[["rsrp1", "smarce1"]].to_numpy()

# sel = VarianceThreshold(threshold=0.005)
# df = sel.fit_transform(df)
# df.shape
# df = SelectKBest(chi2, k=500).fit_transform(df, [1,1,1,1,1,-1,-1,-1,-1,-1])
# df.shape

df = d_53
features = d_53_features
sample = ["b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5"]
pair = { 0:5, 1:6, 2:7, 3:8, 4:9 }
coef_df = pandas.DataFrame()
good_i_53 = []
success = False
for i in range(1, 75):
    # if i % 5 == 0:
    print("i: " + str(i))
    for j in range(1, i+1):
        print("j: " + str(j))
        all_train_score = []
        all_test_score = []
        for b_pair, c_pair in pair.items():
            # print("Predicted Pair: " + sample[b_pair] + " and " + sample[c_pair])
            this_pair = [b_pair, c_pair]
            xtrain = [ x for x in range(0,10) if x not in this_pair ]
            xtrain = df[xtrain]
            xtest = df[this_pair]
            ytrain = numpy.array([1, 1, 1, 1, -1, -1, -1, -1])
            ytest = numpy.array([1, -1])
            rc = LogisticRegression()
            scaler = StandardScaler()
            a = scaler.fit(xtrain)
            xtrain = scaler.transform(xtrain)
            xtest = scaler.transform(xtest)
            # rc = RidgeClassifier(alpha = 1, normalize = False)
            # rc = RandomForestClassifier(n_estimators = 1000, max_depth=100, random_state=0)
            # My Feature Selection
            greater_df = numpy.array([xtrain[x, :] > xtrain[x+4,:] for x in range(0,4)])
            greater_idx = numpy.argwhere(numpy.logical_and(greater_df[0, :], numpy.logical_and(greater_df[1, :], numpy.logical_and(greater_df[2, :], greater_df[3, :]))))
            smaller_df = numpy.array([xtrain[x, :] < xtrain[x+4,:] for x in range(0,4)])
            smaller_idx = numpy.argwhere(numpy.logical_and(smaller_df[0, :], numpy.logical_and(smaller_df[1, :], numpy.logical_and(smaller_df[2, :], smaller_df[3, :]))))
            all_idx = greater_idx.flatten().tolist() + smaller_idx.flatten().tolist()
            xtrain = xtrain[0:8, all_idx]
            xtest = xtest[0:2, all_idx]
            # minmax = MinMaxScaler().fit(xtrain)
            # xtrain = minmax.transform(xtrain)
            # xtest = minmax.transform(xtest)
            # this_features = features[all_idx]
            # print(len(all_idx))
            b_sum = xtrain[0:4, :].sum(axis=0)
            c_sum = xtrain[5:8, :].sum(axis=0)
            dif = abs(b_sum - c_sum)
            dif_idx = (-dif).argsort()
            xtrain = xtrain[0:8, dif_idx[range(0, i)]]
            xtest = xtest[0:2, dif_idx[range(0, i)]]
            # this_features = this_features[dif_idx[range(0, i)]]
            # a = scaler.fit(xtest)  # according to the internet I shouldn't do this step
            a = rc.fit(xtrain, ytrain)
            # importance = rc.feature_importances_
            importance = numpy.abs(rc.coef_)
            sort_idx = (-importance).argsort()
            xtrain = xtrain[0:8, sort_idx[0, range(0,j)]]
            xtest = xtest[0:2, sort_idx[0, range(0,j)]]
            # print(importance[0,sort_idx[0, range(0,5)]])
            # xtrain = xtrain[0:8, sort_idx[0:j]]
            # xtest = xtest[0:2, sort_idx[0:j]]
            # this_features = this_features[sort_idx[0, range(0,j)]]
            a = rc.fit(xtrain, ytrain)
            # rc.predict(xtest)
            # rc.predict_proba(xtest)
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
            # print("Accuracy on test set: ", str(test_score))
            if test_score != 1:
                break
            else:
                all_test_score.append(test_score)
            if len(all_test_score) == 5:
                good_i_53.append([i,j])
                if not success:
                    print("Found a winner!")
                    print(i)
                    print(j)
                    success = True
                    break
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

