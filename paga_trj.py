import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import numpy as np
import anndata2ri
from rpy2.robjects import r

# sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.logging.print_versions()
# results_file = './write/paul15.h5ad'
# sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3))  # low dpi (dots per inch) yields small inline figures

anndata2ri.activate()
adata = r('as("C:/Users/miles/Downloads/d_tooth/data/tj_jaw_sce.RDS", "SingleCellExperiment")')
