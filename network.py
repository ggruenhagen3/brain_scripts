import networkx as nx
import csv
import numpy
import pandas as pd
pd_df = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_j.csv')
# pd_df2 = pd_df.iloc[:,1:]  # if csv has rownames
dj = pd_df.to_numpy()
dj_orig = dj
dj_lower = numpy.tril(dj)
numpy.fill_diagonal(dj_lower, 0)  # remove self to self edges
G = nx.convert_matrix.from_numpy_matrix(dj_lower)  # killed on interactive

new_labels = {}
for x in range(0, pd_df.shape[0]-1):
    new_labels[x] = pd_df.columns[x]

G = nx.relabel.relabel_nodes(G, new_labels)
nx.readwrite.gpickle.write_gpickle(G, "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_j.pickle")

list(G.nodes)[1:5]

import networkx as nx
import csv
import numpy
import pandas as pd
G = nx.readwrite.gpickle.read_gpickle("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_j.pickle")
degrees = sorted(d for n, d in G.degree())
sort_order = numpy.argsort(degrees)
nodes_sorted = numpy.array(G.nodes)[sort_order]

# Sum the weights of the edges into/out of each node
gene_weights = {}
for n, nbrs in G.adj.items():
    sum_w = 0
    for nbr, eattr in nbrs.items():
        sum_w += eattr['weight']
    gene_weights[n] = sum_w

# Print the top 5
gene_weights_sort = dict(sorted(gene_weights.items(), key=lambda item: item[1]))
gene_weights_sort = list(gene_weights_sort.keys())
gene_weights_sort.reverse()
gene_weights_sort[0:5]

# Add the pagerank to the node data
G = nx.readwrite.gpickle.read_gpickle("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_j.pickle")
pr = nx.pagerank(G)
for key, value in pr.items():
    G.nodes[key]["pr"] = value
nx.readwrite.gpickle.write_gpickle(G, "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_c_j_w_pr.pickle")

pr_b = nx.readwrite.gpickle.read_gpickle("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_b_pr.pickle")
pr_c = nx.readwrite.gpickle.read_gpickle("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_c_pr.pickle")
# for key, value in pr.items():
#     G.nodes[key]["pr"] = value
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/pr_b.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in pr_b.items():
       writer.writerow([key, value])
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/pr_c.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in pr_c.items():
       writer.writerow([key, value])



pos = nx.nx_agraph.graphviz_layout(G)

nx.draw(G, pos=pos)
write_dot(G, 'file.dot')

"""
RNA Velocity
"""
import scvelo as scv
from matplotlib import pyplot as plt
scv.set_figure_params()
adata = scv.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-B1/velocyto/JTS07-B1.loom", cache=True)
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

import loompy
loompy.combine(["/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-B1/velocyto/JTS07-B1.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-B2/velocyto/JTS07-B2.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-B3/velocyto/JTS07-B3.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-B4/velocyto/JTS07-B4.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-B5/velocyto/JTS07-B5.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-C1/velocyto/JTS07-C1.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-C2/velocyto/JTS07-C2.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-C3/velocyto/JTS07-C3.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-C4/velocyto/JTS07-C4.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-C5/velocyto/JTS07-C5.loom"],
               "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/all_velo/all_velo.loom")
adata = scv.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/all_velo/all_velo.loom", cache=True)
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
seurat_data = scv.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_seurat.loom")

scv.tl.umap(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', save = "veloplot.svg")

seurat_data = scv.read("b1.loom")
merged = scv.utils.merge(adata, seurat_data)
scv.pl.velocity_embedding_stream(merged, basis='umap_cell_embeddings', color='seurat_clusters', save = "_all_merged.svg")
scv.pl.velocity_embedding_grid(merged, basis='umap_cell_embeddings', color='seurat_clusters', save = "_all_merged_grid.svg")
scv.pl.velocity_embedding(merged, basis='umap_cell_embeddings', color='seurat_clusters', save = "_all_merged_cell.svg")

ad_sub=merged[merged.obs['seuratclusters15']=='0',:]
scv.pl.velocity_embedding_stream(ad_sub, basis='umap_cell_embeddings', color='seurat_clusters', save = "_all_merged0_stream.svg")
scv.pl.velocity_embedding_grid(ad_sub, basis='umap_cell_embeddings', color='seurat_clusters', save = "_all_merged0_grid.svg")

scv.pp.filter_and_normalize(ad_sub)
scv.pp.moments(ad_sub)
scv.tl.velocity(ad_sub, mode='stochastic')
scv.tl.velocity_graph(ad_sub)

merged = scanpy.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/all_velo/all_merged_w_seurat_dynamical.h5ad")
scv.tl.latent_time(merged)
merged.obs['root_cells'] = 1 - merged.obs['root_cells']
merged.obs['end_points'] = 1 - merged.obs['end_points']
# merged.obs.loc[merged.obs['seuratclusters53'] != '12']['root_cells'] = 0
# merged.obs.loc[merged.obs['seuratclusters53'] == '12']['root_cells'] = 1
scv.tl.recover_latent_time(merged)
ad_sub=merged[merged.obs['seuratclusters15']=='0',:]
# scv.tl.recover_latent_time(ad_sub)
scv.pl.scatter(ad_sub, color='latent_time', color_map='Spectral_r', size=80, basis = 'umap_cell_embeddings', save = '_dynamical_time0_flip.svg')

import scanpy
scanpy.AnnData.write_loom(merged, "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/all_velo/all_merged_w_seurat.loom")
'write/pbmc3k.h5ad'
merged.write("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/all_velo/all_merged_w_seurat.h5ad")

nodes = list(G.nodes)
for node in nodes[0:(len(nodes)-10)]:
    G.remove_node(node)
nx.readwrite.gpickle.write_gpickle(G, "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/test.pickle")

pos = nx.nx_agraph.graphviz_layout(G)

"""
Tooth RNA Velocity
"""
import loompy
import scvelo as scv
loompy.combine(["/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/bs/MZ/UJ_bcl/velocyto/UJ_bcl.loom",
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/bs/MZ/LJ_bcl/velocyto/LJ_bcl.loom"],
                "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/bs/MZ/all_jaw.loom")
adata = scv.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/bs/MZ/all_jaw.loom", cache=True)
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
seurat_data = scv.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/jaw_loom.loom")

merged = scv.utils.merge(adata, seurat_data)
scv.pl.velocity_embedding_stream(merged, basis='umap_cell_embeddings', color='seurat_clusters', save = "jaw_veloplot.svg")
scv.pl.velocity_embedding_grid(merged, basis='umap_cell_embeddings', color='seurat_clusters', save = "jaw_veloplot_grid.svg", arrow_size = 2)
scv.pl.velocity_embedding(merged, basis='umap_cell_embeddings', color='seurat_clusters', save = "jaw_veloplot_cell.svg", arrow_size = 2)

