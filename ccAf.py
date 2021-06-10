# conda activate ccAf2
import scanpy as sc
import ccAF
import pandas

pd_df = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/rownames_to_human_to_humanID.csv")

# TJ_JAW
# Convert names
obj = sc.read_loom('/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/tj_jaw.loom')
obj.var_names = pd_df["hgnc.id"].fillna("").tolist()

# Predict cell cycle phase labels
predictedLabels = ccAF.ccAF.predict_labels(obj)
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/tj_jaw_ccAf.txt', 'w') as filehandle:
    filehandle.writelines("%s\n" % place for place in predictedLabels)

# Jaw
obj = sc.read_loom('/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/jaw_loom.loom')
obj.var_names = pd_df["hgnc.id"].fillna("").tolist()

# Predict cell cycle phase labels
predictedLabels = ccAF.ccAF.predict_labels(obj)
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/results/jaw_ccAf.txt', 'w') as filehandle:
    filehandle.writelines("%s\n" % place for place in predictedLabels)
