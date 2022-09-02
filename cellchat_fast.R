# Read Input ===================================================================
# this.run = 1; do.down = T; is.real = F; num.perms = 100;
args = commandArgs(trailingOnly=TRUE)
this.run  = as.numeric(args[1])
do.down   = as.logical(args[2])
is.real   = as.logical(args[3])
num.perms = as.numeric(args[4])
set.seed(this.run)
message(paste0("Initializng run with parameters: this.run=", this.run, ", do.down=", do.down, ", num.perms=", num.perms, "."))

# Load Libraries ===============================================================
suppressMessages(library('CellChat',  quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('patchwork', quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('stringr',   quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('parallel',  quietly = T, warn.conflicts = F, verbose = F))
options(stringsAsFactors = FALSE)

# Load Data ====================================================================
message("Loading bb...")
setwd("~/scratch/brain/cellchat/")
gene_info = read.table("gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
combined = readRDS("bb_demux_012422.rds")
labels = read.csv("primary_w_rgc_labels.csv")
combined$label = labels$x
meta = data.frame(label = combined$label, row.names = colnames(combined))
message("Done.")

# Downsample ===================================================================
if (do.down) {
  message("Doing Downsampling...")
  down_sample_fnc = function() { unlist(lapply(unique(combined$label), function(x) sample(colnames(combined)[which(combined$label == x)], 49) )) }
  down_cell_list  = mclapply(1:num.perms, function(x) down_sample_fnc(), mc.cores = 24)
  message("Done.")
} else { message("Not Downsampling.") }

# Permutations =================================================================
getLabels = function(x) {
  this.labels = c()
  if (is.real) {
    if (do.down) { this.labels = combined$label[down_cell_list[[x]]]         } else { this.labels = combined$label         }
  } else {
    if (do.down) { this.labels = sample(combined$label[down_cell_list[[x]]]) } else { this.labels = sample(combined$label) }
  }
  return (this.labels)
}

if (is.real) { message("Gather Cell Labels...") } else { message("Permuting Data...") }
label_list = mclapply(1:num.perms, function(x) getLabels(x), mc.cores = 24)
message("Done.")

# Human Object =================================================================
message("Creating a Human Object...")
mz.df = data.frame(mz = rownames(combined), human = gene_info$human[match(rownames(combined), gene_info$mzebra)])
mz.df$rowsums = rowSums(combined@assays$RNA@data)
mz.df = mz.df[order(-mz.df$rowsums),]
mz.df = mz.df[which(mz.df$rowsums != 0),]

mz.df = mz.df[!duplicated(mz.df$human),]
data.input = as.matrix(combined@assays$RNA@data[mz.df$mz,])
rownames(data.input) = mz.df$human
message("Done.")

rm(combined) # delete original Seurat object to save memory

# Cell Chat ====================================================================
CellChatWeights = function(x) {
  if (do.down) { this.cells = down_cell_list[[x]] } else { this.cells = colnames(data.input) }
  this.meta = data.frame(label = label_list[[x]], row.names = this.cells)
  
  cellchat = createCellChat(object = data.input[,this.cells], meta = this.meta, group.by = "label")
  cellchat = addMeta(cellchat, meta = this.meta)
  cellchat = setIdent(cellchat, ident.use = "label")
  
  cellchat@DB = CellChatDB.human
  cellchat = subsetData(cellchat)
  # invisible(capture.output( cellchat = identifyOverExpressedGenes(cellchat)        ))
  # invisible(capture.output( cellchat = identifyOverExpressedInteractions(cellchat) ))
  cellchat = identifyOverExpressedGenes(cellchat)
  cellchat = identifyOverExpressedInteractions(cellchat)
  
  cellchat = computeCommunProb(cellchat, type =  "truncatedMean", trim = 0.1, population.size = F)
  df.net_lig_recept = subsetCommunication(cellchat) 
  cellchat = aggregateNet(cellchat)
  cellchat = filterCommunication(cellchat, min.cells = 10)
  net_weight = data.frame(cellchat@net$weight)
  net_weight_vect = unlist(net_weight)
  name_rep = rep(rownames(net_weight), ncol(net_weight))
  names(net_weight_vect) = paste0(name_rep, ".", sort(name_rep))
  return(net_weight_vect)
}

message("Running cellchat (this while take awhile)...")
if (do.down) { num.parallel.jobs = 4 } else { num.parallel.jobs = 2 }
# onerun = suppressMessages(CellChatWeights(1))
sink(file="~/scratch/brain/cellchat_sink.txt")
run_outs = mclapply(1:num.perms, function(x) suppressMessages(CellChatWeights(x)), mc.cores = num.parallel.jobs)
sink()

out = do.call('cbind', run_outs)
colnames(out) = c("run", 1:num.perms)
out[, c("clust1", "clust2")] = reshape2::colsplit(names(run_outs[[1]]), "\\.", c("1", "2"))
out = out[, c(num.perms+1, num.perms+2, 1:num.perms)]
message("Done.")

# Save Output ==================================================================
message("Writing Output...")
todays.date = stringr::str_split(Sys.Date(), pattern = "-")[[1]]
todays.date = paste0(todays.date[2], todays.date[3], substr(todays.date[1], 3, 4))
out.str = paste0("~/scratch/brain/results/cellchat/primary/cellchat_", ifelse(do.down, "downsampled_", "full_"), ifelse(is.real, "real_", "perm_"), num.perms, "nruns_run", this.run, ".csv")
write.csv(out, out.str)
message("Done.")
message("All Done.")
