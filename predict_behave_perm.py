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
import multiprocessing

def permSingleRun(i):
    # Set Constants for Number of Features to Select
    num_cum_dif = 117
    num_top_model = 4

    # Read Data and Separate Train from Test
    pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/pseudo_subsample_ml/pseudo_" + str(i) + ".csv", index_col=0)
    train_names = [x for x in list(pd_df.index) if "train" in x]  # the subsamples that are in the test data
    test_names = [x for x in list(pd_df.index) if "test" in x]  # the subsamples that are in the test data
    xtrain = pd_df.loc[train_names]
    xtest = pd_df.loc[test_names]

    # Find Genes with the largest cummulative difference between bhve and ctrl
    bhve_names = [x for x in list(pd_df.index) if "b" in x]
    ctrl_names = [x for x in list(pd_df.index) if "c" in x]
    b_sum = xtrain.loc[bhve_names].sum(axis=0)
    c_sum = xtrain.loc[ctrl_names].sum(axis=0)
    dif = abs(b_sum - c_sum)
    dif_idx = (-dif).argsort()
    xtrain = xtrain.iloc[:, dif_idx[0:num_cum_dif]]
    xtest = xtest.iloc[:, dif_idx[0:num_cum_dif]]

    # Select the Top Features to the Model
    rc = LogisticRegression(C=1)
    a = rc.fit(xtrain, ytrain)
    importance = numpy.abs(rc.coef_)
    sort_idx = (-importance).argsort()  # argsort give the indexes which would sort the array
    xtrain = xtrain.iloc[:, sort_idx[0, 0:num_top_model]]
    xtest = xtest.iloc[:, sort_idx[0, 0:num_top_model]]

    # Refit the Model and Predict
    a = rc.fit(xtrain, ytrain)
    test_score = rc.score(xtest, ytest)

    return(test_score)


# 100 Permutations
with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
    pool_res = pool.map(permSingleRun, range(0, 200))
    perm_df = pandas.DataFrame(pool_res, columns = list(range(0, 200)))


