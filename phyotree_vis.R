library(tidyverse)
library(ape)
library(ggmap)
library(phyloseq)

my_path = "C:/Users/miles/Downloads/"

# tree <- ape::read.tree(paste0(my_path, "pit_test.newick"))
# tree <- ape::read.tree(paste0(my_path, "pit_final_nj.newick"))
# tree <- ape::read.tree(paste0(my_path, "castle_rock_pit_nj.newick"))
# tree <- ape::read.tree(paste0(my_path, "castle_rock_pit_JTS09_nj.newick"))
# tree <- ape::read.tree(paste0(my_path, "castle_rock_pit_JTS09_50_nj.newick"))
# tree <- ape::read.tree(paste0(my_path, "snphylo.output.ml.tree"))
# tree <- ape::read.tree(paste0(my_path, "snphylo.output2.ml.tree"))
# tree <- ape::read.tree(paste0(my_path, "snphylo.output.no.JTS09.ml.tree"))
tree <- ape::read.tree(paste0(my_path, "snphylo.pvc.ml.tree"))

# Optionally set an outgroup.
# tree <- root(tree,outgroup = "outgroup", resolve.root = T)

treeSegs <- phyloseq::tree_layout(
  phyloseq::phy_tree(tree),
  ladderize = T
)

treeSegs$edgeDT <- treeSegs$edgeDT  %>%
  dplyr::mutate(edge.length =
                  ifelse(edge.length < 0, 0, edge.length)
                , xright = xleft + edge.length
  )
edgeMap = aes(x = xleft, xend = xright, y = y, yend = y)
vertMap = aes(x = x, xend = x, y = vmin, yend = vmax)
labelMap <- aes(x = xright+0.0001, y = y, label = OTU)

# ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) + 
#   geom_segment(vertMap, data = treeSegs$vertDT) +
#   geom_text(labelMap, data = dplyr::filter(treeSegs$edgeDT, !is.na(OTU)), na.rm = TRUE, hjust = -0.05) +
#   ggmap::theme_nothing() + 
#   scale_x_continuous(limits = c(
#     min(treeSegs$edgeDT$xleft)-0.15,
#     max(treeSegs$edgeDT$xright)+0.15
#   ),
#   expand = c(0,0))


pit_species = c("CV", "NP", "FR", "TP", "DK", "DC", "TI", "MS", "ML", "AB", "AC")
castle_species = c("TF", "CL", "CM", "MA", "CN", "OA", "MC", "NO", "TC")
rock_species = c("MZ", "LF", "CA", "CO", "MP", "PC", "GM")
JTS09_species = c("1B11", "1B18", "1B4", "1B5", "1C11", "1C17", "1C4", "1C5", "2B10", "2B13", "2B17", "2B19", "2C10", "2C14", "2C18", "2C19", "3B2", "3B6", "3B7", "3B9", "3C2", "3C6", "3C7", "3C9", "4B12", "4B14", "4B25", "4C12", "4C13", "4C25", "5B22", "5B23", "5B24", "5B26", "5C22", "5C23", "5C24", "5C26")
labelMap <- aes(x = xright+0.0001, y = y, label = OTU, color = bower, fill = bower)
treeSegs$edgeDT$OTU = str_replace(treeSegs$edgeDT$OTU, "_all", "")
treeSegs$edgeDT$bower = treeSegs$edgeDT$OTU
treeSegs$edgeDT$bower[which(treeSegs$edgeDT$OTU %in% pit_species)] = "pit"
treeSegs$edgeDT$bower[which(treeSegs$edgeDT$OTU %in% castle_species)] = "castle"
treeSegs$edgeDT$bower[which(treeSegs$edgeDT$OTU %in% JTS09_species)] = "JTS09"
treeSegs$edgeDT$bower[which(treeSegs$edgeDT$OTU %in% rock_species)] = "rock"
ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) + 
  geom_segment(vertMap, data = treeSegs$vertDT) +
  geom_label(labelMap, data = dplyr::filter(treeSegs$edgeDT, !is.na(OTU)), na.rm = TRUE, hjust = -0.05, alpha = 0.2) +
  # geom_text(aes(x = (xright + xleft) / 2, y = y, label = edge.length, vjust = -0.2), data = dplyr::filter(treeSegs$edgeDT, !is.na(OTU)), na.rm = TRUE) +
  ggmap::theme_nothing() + 
  scale_x_continuous(limits = c(
    min(treeSegs$edgeDT$xleft)-0.15,
    max(treeSegs$edgeDT$xright)+0.15
  ),
  expand = c(0,0))
