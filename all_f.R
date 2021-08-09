library("stringr")
library("ggplot2")
library("biomaRt")
library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("RColorBrewer")
library("dplyr")
library("gplots")
library("viridis")
library("reshape2")
library("edgeR")
library("data.table")
library("tidyr")
library("BSDA")
library("httr")
library("ggtext")
library("colourvalues")
library("colorspace")
library("ggrepel")
# library("ggforce")
httr::set_config(config(ssl_verifypeer = FALSE))

####################
# Helper Functions #
####################

clipboard <- function(x, sep="\t", row.names=FALSE, col.names=TRUE){
  con <- pipe("xclip -selection clipboard -i", open="w")
  write.table(x, con, sep=sep, row.names=row.names, col.names=col.names, quote = F)
  close(con)
}

pct_dif_avg_logFC = function(obj, cells.1, cells.2) {
  #' Find the Average Log FC and Percent Difference for Every Gene
  #' @param obj Seurat object
  #' @param cells.1 vector of first group of cells
  #' @param cells.2 vector of second group of cells
  both_cells = unique(c(cells.1, cells.2))
  non_zero_genes = rownames(obj)[which(rowSums(obj@assays$RNA@counts[,both_cells]) > 0)]
  
  # Find Percent of Gene+ in both groups of cells
  mat1 = obj@assays$RNA@counts[non_zero_genes, cells.1]
  mat1[which(mat1 > 1)] = 1
  num_pos_cells1 = rowSums(mat1)
  mat2 = obj@assays$RNA@counts[non_zero_genes, cells.2]
  mat2[which(mat2 > 1)] = 1
  num_pos_cells2 = rowSums(mat2)
  pct_pos_cells1 = (num_pos_cells1/length(cells.1)) * 100
  pct_pos_cells2 = (num_pos_cells2/length(cells.2)) * 100
  pct_dif = pct_pos_cells1 - pct_pos_cells2
  
  # Find Average Log FC
  # Seurat
  mean_exp1 = rowMeans(expm1(obj@assays$RNA@data[non_zero_genes, cells.1]))
  mean_exp2 = rowMeans(expm1(obj@assays$RNA@data[non_zero_genes, cells.2]))
  avg_logFC = log(mean_exp1 + 1) - log(mean_exp2 + 1)
  
  # # Zack 
  # mean_exp1 = rowSums(bb@assays$RNA@counts[non_zero_genes,cells.1]) / sum(bb$nCount_RNA[cells.1])
  # mean_exp2 = expm1(rowSums(bb@assays$RNA@counts[non_zero_genes,cells.2]) / sum(bb$nCount_RNA[cells.2]))
  # avg_logFC = log(mean_exp1 + 1) - log(mean_exp2 + 1)
  
  df = data.frame(genes = non_zero_genes, avg_logFC = avg_logFC, pct.1 = pct_pos_cells1, pct.2 = pct_pos_cells2, pct_dif = pct_dif, num.1 = num_pos_cells1, num.2 = num_pos_cells2)
  
  return(df)
}

r_to_p = function(r1, r2, n1, n2) {
  # Compare Two Correlation Values using Fisher's Z Transformation Method.
  #' @param r1 correlation 1
  #' @param r2 correlation 2
  #' @param n1 number of samples used for correlation 1
  #' @param n2 number of samples used for correlation 2
  z1 = .5 * (log(1+r1) - log(1-r1))
  z2 = .5 * (log(1+r2) - log(1-r2))
  z3 =(z1 - z2) / sqrt( (1 / (n1 - 3)) + (1 / (n2 - 3)) )
  p = 2*pnorm(-abs(z3))
  return(p)
}

changeClusterID = function(clust_vect, clust_level = NULL, returnFactor = F) {
  #' Convert old cluster labels to the new labels
  #' @param clust_vect vector of old cluster labels
  #' @param clut_level 15 or 53 cluster level? If NA, detect automatically
  #' @param returnFactor return vector as a factor?
  #' @return new_vect a vector of new cluster labels
  
  # If not supplied, determine the cluster level of the input
  clust_vect = as.numeric(as.vector(clust_vect))
  if (is.null(clust_level)) {
    max_clust = max(clust_vect)
    clust_level = 15
    if (max_clust > 14)
      clust_level = 53
    print(paste0("No cluster level input. Assuming ", clust_level, " level because the maximum cluster number found in the input was ", max_clust, "."))
  }
  
  # Build Both 15 and 53 level converters
  convert15 = data.frame(old = 0:14, new = c("8_Ex", "9_Ex", "4_In", "15_In/Ex", "1_Astro/MG", "10_Ex", "5_In", "11_Ex", "6_In", "2_OPC/Oligo", "12_Ex", "13_Ex", "14_Ex", "3_Peri", "7_In"))
  convert53 = data.frame(old = 0:52, new = c("4.1_In", "10.1_Ex", "15.1_In/Ex", "9.1_Ex", "8.1_Ex", "1.1_Astro", "6_In", "5.1_In", "9.2_Ex", "8.2_Ex", "15.2_In", "11.1_Ex", "8.3_Ex", "8.4_Ex", "9.3_Ex", "4.2_In", "8.5_Ex", "5.2_In", "8.6_Ex", "8.7_Ex", "1.2_Astro", "4.3_In", "4.4_In", "9.4_Ex", "9.5_Ex", "8.8_Ex", "9.6_Ex", "4.5_In", "12_Ex", "8.9_Ex", "10.2_Ex", "2.1_OPC", "15.3_In", "11.2_Ex", "15.4_In", "4.6_In", "9.7_Ex", "13_Ex", "14_Ex", "4.7_In", "11.3_Ex", "9.8_Ex", "8-9_Ex", "15.5_In/Ex", "4.8_In", "1.3_MG", "2.2_Oligo", "15.6_Ex", "8.10_Ex", "8.11_Ex", "3_Peri", "15.7_Ex", "7_In"))
  
  # Changed In and Ex to GABA and GLUT
  convert15$new = gsub("In", "GABA", convert15$new)
  convert15$new = gsub("Ex", "Glut", convert15$new)
  convert53$new = gsub("In", "GABA", convert53$new)
  convert53$new = gsub("Ex", "Glut", convert53$new)
  
  # Allow converters to sort numerically
  convert15 = cbind(convert15, colsplit(convert15$new, pattern = "_", names = c("new.num", "new.gaba")))
  convert15 = convert15[order(convert15$new.num, decreasing = F),]
  convert53 = cbind(convert53, colsplit(convert53$new, pattern = "_", names = c("new.num", "new.gaba")))
  convert53 = cbind(convert53, colsplit(convert53$new.num, pattern = "\\.", names = c("new.num1", "new.num2")))
  convert53$new.num1[which(convert53$new == "8-9_Glut")] = 8
  convert53$new.num2[which(convert53$new == "8-9_Glut")] = 99
  convert53 = convert53[order(as.numeric(convert53$new.num1), as.numeric(convert53$new.num2), decreasing = F),]
  
  # Select correct converter
  converter = convert53
  if (clust_level == 15)
    converter = convert15
  
  # Convert Data
  new_vect = as.vector(converter$new[match(clust_vect, converter$old)])
  
  if (returnFactor)
    new_vect = factor(new_vect, levels = converter$new)
  
  print("Conversion Successful.")
  return(new_vect)
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                          draw_group = function(self, data, ..., draw_quantiles = NULL) {
                            data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                            grp <- data[1, "group"]
                            newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                            newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                            newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                            
                            if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                              stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                        1))
                              quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                              aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                              aesthetics$alpha <- rep(1, nrow(quantiles))
                              both <- cbind(quantiles, aesthetics)
                              quantile_grob <- GeomPath$draw_panel(both, ...)
                              ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                            }
                            else {
                              ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                            }
                          })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

bvcVis = function(obj, feature, myslot = "data", mode = "box", meta = "sample", only.pos = F, cells.use = NULL, cell.alpha = 0.08, add_title = "") {
  #' Behave vs Control Replicate Boxplot Visualization
  #' 
  #' @param obj Seurat object
  #' @param cell.alpha alpha value for the cells (the geom jitter points)
  #' @param feature feature to visualize
  #' @param myslot slot to pull data from
  #' @param mode mode of visualization ("box" for boxlpot, "violin" for violinplot, "violin_split" for violin plot split by meta)
  #' @param meta meta variable to split data by
  #' @param only.pos plot only the positive cells
  #' @param cells.use the cells to plot
  #' @param add_title input string to add to the title
  
  if (feature %in% rownames(obj)) {
    exp = GetAssayData(obj, assay = "RNA", slot=myslot)
    values = exp[feature,]
  } else {
    values = as.numeric(as.vector(obj@meta.data[c(feature)][,c(feature)]))
  }
  df = data.frame(values, sample = obj@meta.data[c(meta)][,c(meta)]) 
  rownames(df) = colnames(obj)
  
  if (! is.null(cells.use))
    df = df[cells.use,]
  if (only.pos)
    df = df[which(df$values > 0),]
  
  if (length(unique(obj@meta.data[c(meta)][,c(meta)])) == 2)
    col_pal = c("#F2444A","#0077b6")
  else
    col_pal = c("#9d0208", "#d00000", "#dc2f02", "#e85d04", "#f48c06", "#03045e", "#023e8a", "#0077b6", "#0096c7", "#00b4d8")
  # col_pal = c("#9d0208", "#d00000", "#dc2f02", "#e85d04", "#f48c06", "#7400b8", "#6930c3", "#5e60ce", "#5390d9", "#4ea8de")
  # col_pal = c("#F2444A","#F25C44", "#F27A44", "#F2B644", "#F2D444", "#023e8a", "#0077b6", "#0096c7", "#00b4d8", "#48cae4")
  # col_pal = c("#F2444A","#F25C44", "#F27A44", "#F2B644", "#F2D444", "#023e8a", "#0077b6", "#0096c7", "#00b4d8", "#48cae4")
  # col_pal = c("#F2444A","#F25C44", "#F27A44", "#F2B644", "#F2D444", "#22577a", "#38a3a5", "#57cc99", "#80ed99", "#918ef4")
  # col_pal = c("#F2444A","#F25C44", "#F27A44", "#F2B644", "#F2D444", "#7400b8", "#6930c3", "#5e60ce", "#5390d9", "#4ea8de")
  my_title = paste(feature, add_title)
  
  if (mode == "box")
    p = ggplot(df, aes(x=sample, y = values, fill=sample, color=sample)) + geom_boxplot(alpha=0.6) + geom_jitter(position=position_dodge2(width = 0.7), alpha = cell.alpha) + scale_color_manual(values=col_pal) + scale_fill_manual(values=col_pal) + ylab("Normalized Expression") + xlab("") + ggtitle(my_title) + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
  if (mode == "violin")
    p = ggplot(df, aes(x=sample, y = values, fill=sample, color=sample)) + geom_violin(alpha=0.6) + geom_jitter(position=position_dodge2(width = 0.6), alpha = cell.alpha) + scale_color_manual(values=col_pal) + scale_fill_manual(values=col_pal) + ylab("Normalized Expression") + xlab("") + ggtitle(my_title) + theme(plot.title = element_text(hjust = 0.5)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="#495057", fun.args = list(mult = 1)) + NoLegend()
  if (mode == "violin_split") {
    df$pair = substr(as.character(as.vector(df$sample)), 2, 2)
    if (meta == "cond")
      p = ggplot(df, aes(x=1, y = values, fill=sample, color=sample)) + geom_split_violin(alpha=0.6) + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)) + scale_color_manual(values=col_pal, name = "Condition") + scale_fill_manual(values=col_pal, name = "Condition") + ylab("Normalized Expression") + xlab("") + ggtitle(my_title) + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_blank(), axis.ticks=element_blank())
    else
      p = ggplot(df, aes(x=pair, y = values, fill=sample, color=sample)) + geom_split_violin(alpha=0.6) + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)) + scale_color_manual(values=col_pal, name = "Sample") + scale_fill_manual(values=col_pal, name = "Sample") + ylab("Normalized Expression") + xlab("Sample Pair") + ggtitle(my_title) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(p)
}

myDensityPlot = function(obj, feature1, feature2, cells.use = NULL, myslot = "data", alpha_vect = NULL, my.split.by = NULL, my.pt.size = 1) {
  #' Make a FeaturePlot using ggplot. Having this gives me more control over the plot.
  #' 
  #' @param obj Seurat object
  #' @param feature feature to plot
  #' @param myslot slot to pull data from
  #' @param alpha_vect vector of values to use as an alpha parameter (values in same order as cells)
  #' @return ggplot object
  exp = GetAssayData(obj, assay = "RNA", slot="counts")
  values = exp[feature1,] > 0
  values2 = exp[feature2,] > 0
  
  df = as.data.frame(obj@reductions$umap@cell.embeddings)
  df$pos = 0
  df$pos[which(values)] = 1
  df$pos[which(values2)] = 2
  df$pos[which(values & values2)] = 3
  df$ident = Idents(obj)
  df$alpha = 0.5
  df$alpha[which(df$pos == feature1)] = 0.8
  
  hull <- df %>% group_by(pos) %>%
    slice(c(chull(UMAP_1,UMAP_2),
            chull(UMAP_1,UMAP_2)[1]))
  
  df = df[order(df$pos),]
  df$pos = plyr::revalue(as.character(df$pos), replace = c("1" = feature1, "2" = feature2, "3" = paste(feature1, "+", feature2), "0" = "none"))
  # p = ggplot(df, aes(UMAP_1, UMAP_2, col = pos)) + geom_point(size = my.pt.size) + stat_ellipse(level = 0.95) + theme_classic() + theme(plot.title = element_text(hjust = 0.5))
  p = ggplot(data=df, mapping = aes(UMAP_1, UMAP_2, col = pos)) + geom_point(size = my.pt.size, alpha = 1)  + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = c("#2a9d8f", "#e76f51", "#264653", "lightgray"))
  # p = p + geom_density_2d(data=df[which( df$pos == paste(feature1, "+", feature2) ),], mapping=aes(UMAP_1, UMAP_2, col = pos), size = 1, contour_var = "density", position = "identity", bins = 5, alpha = 0.8)
  # p = p + geom_density_2d_filled(data=df[which( df$pos == paste(feature1, "+", feature2) ),], mapping=aes(UMAP_1, UMAP_2, col = pos), contour_var = "ndensity", position="identity", bins = 5, alpha = 0.4)
  
  # I haven't gotten this to work, but this will let me use two different color scales
  # for the two different genes
  # https://stackoverflow.com/questions/20474465/using-different-scales-as-fill-based-on-factor
  # # plot with red fill
  # p = ggplot(data=df, mapping = aes(UMAP_1, UMAP_2, col = pos)) + geom_point(size = my.pt.size, alpha = 1)  + theme_classic() + theme(plot.title = element_text(hjust = 0.5))
  # p1 <- p + stat_density2d(data = df[which(df$pos == feature1 | df$pos == feature2),], aes(UMAP_1, UMAP_2, color = as.factor(pos), fill = ..level..), alpha = 0.3, geom = "polygon") +
  #   scale_fill_continuous(low = "grey", high = "red", space = "Lab", name = "g = 0") +
  #   scale_colour_discrete(guide = FALSE) +
  #   theme_classic() + ggtitle("p1")
  # 
  # # plot with blue fill
  # p2 <- p + stat_density2d(data = df[which(df$pos == feature1 | df$pos == feature2),], aes(UMAP_1, UMAP_2, color = as.factor(pos), fill = ..level..), alpha = 0.3, geom = "polygon") +
  #   scale_fill_continuous(low = "grey", high = "blue", space = "Lab", name = "g = 0") +
  #   scale_colour_discrete(guide = FALSE) +
  #   theme_classic() + ggtitle("p2")
  # 
  # # grab plot data
  # pp1 <- ggplot_build(p1)
  # pp2 <- ggplot_build(p2)$data[[1]]
  # 
  # # replace red fill colours in pp1 with blue colours from pp2 when group is 2
  # pp1$data[[1]]$fill[grep(pattern = "^2", pp2$group)] <- pp2$fill[grep(pattern = "^2", pp2$group)]
  # 
  # # build plot grobs
  # grob1 <- ggplot_gtable(pp1)
  # grob2 <- ggplotGrob(p2)
  # 
  # # build legend grobs
  # leg1 <- gtable_filter(grob1, "guide-box")
  # leg2 <- gtable_filter(grob2, "guide-box")
  # leg <- gtable:::rbind_gtable(leg1[["grobs"]][[1]],  leg2[["grobs"]][[1]], "first")
  # 
  # # replace legend in 'red' plot
  # grob1$grobs[grob1$layout$name == "guide-box"][[1]] <- leg
  # 
  # grid.newpage()
  # grid.draw(grob1)
  # print(hist(c(1,2,3)))
  # print(p1)
  # print(p2)
  # p = "hi"
  
  return(p)
}

multiFeaturePlot2 = function(obj, features, my.pt.size = 1, my.alpha = 0.8, cells.use = NULL, cols = NULL) {
  #' Create a FeaturePlot using multiple markers with different colors
  #' This one is the preferred function, it uses a color blending function, 
  #' but multiFeaturePlot is also acceptable and its method is to overlap points.
  #' @param obj Seurat object
  #' @param features the vector of features to make a FeaturePlot of
  #' @param my.pt.size size of points on plot (optional)
  #' @param my.alpha alpha for all points on plot (optional)
  #' @param cells.use cells to plot (optional)
  #' @param cols vector of colors corresponding to the genes (optional)
  #' @return p UMAP FeaturePlot
  
  # Color Scheme
  lightgray_hex = "#D3D3D3FF"
  if (is.null(cols))
    cols = brewer.pal(9, "Set1")
  
  # Make sure the input features are in the data and get the value of 
  df = as.data.frame(obj@reductions$umap@cell.embeddings)
  df$ident = Idents(obj)
  if (length(features) > 1)
    df$order = colSums(obj@assays$RNA@counts[features,])
  else
    df$order = obj@assays$RNA@counts[features,]
  n = ncol(obj)
  
  i = 1
  for (feature in features) {
    if (feature %in% rownames(obj)) {
      this_col = cols[i]
      if (i == 1) {
        df$feature = obj@assays$RNA@data[feature,]
        m <- grDevices::colorRamp(c(lightgray_hex, this_col))( (1:n)/n )
        df$col = colour_values(df$feature, palette = m)
      } else {
        this_values = obj@assays$RNA@data[feature,]
        m <- grDevices::colorRamp(c("lightgray", this_col))( (1:n)/n )
        this_cols = colour_values(this_values, palette = m)
        names(this_cols) = colnames(obj)
        
        feature_pos = colnames(obj)[which(obj@assays$RNA@counts[feature,] > 0 )]
        newRow = df[feature_pos,]
        newRow$feature = this_values[feature_pos]
        newRow$col = this_cols[feature_pos]
        
        my_mix = function(x)  {
          if (df[x,"col"] == lightgray_hex)
            this_cols[x]
          else
            rgb(mixcolor(alpha = 0.5, hex2RGB(df[x, "col"]), hex2RGB(this_cols[x]))@coords)
        }
        new_cols = unlist(lapply(feature_pos, my_mix))
        df[feature_pos, "col"] = new_cols
      }
      i = i + 1
    } else {
      print("Feature not found in rownames of object nor meta data.")
    }
  }
  
  if (! is.null(cells.use) )
    df = df[cells.use,]
  
  df = df[order(df$order, decreasing = F),]
  print(head(df))
  p = ggplot(df, aes(UMAP_1, UMAP_2)) + geom_point(size = my.pt.size, colour = df$col, alpha = my.alpha) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
  p = LabelClusters(p, "ident", unique(df$ident), repel = F, colour = "black")
  
  # Color title based on genes
  title_str = "<span>"
  for (i in 1:length(features)) {
    if (i == length(features)) {
      title_str = paste0(title_str, "<span style='color:", as.character(cols[i]), ";'>", features[i], "</span>")
    } else {
      title_str = paste0(title_str, "<span style='color:", as.character(cols[i]), ";'>", features[i], "</span> and ")
    }
    if (i %% 3 == 0) {
      title_str = paste0(title_str, "<br></br>")
    }
  }
  title_str = paste0(title_str, "</span>")
  p = p + labs(title = title_str) + theme( plot.title = element_markdown(lineheight = 1.1 ))
  
  return(p)
}

multiFeaturePlot = function(obj, features, my.pt.size = 1) {
  #' Experimental Function to create a FeaturePlot using multiple markers
  #' @param obj Seurat object
  #' @param features the vector of features to make a FeaturePlot of
  
  # cols = rainbow(length(features))
  cols = brewer.pal(9, "Set1")
  # Make sure the input features are in the data and get the value of 
  df = as.data.frame(obj@reductions$umap@cell.embeddings)
  df$ident = Idents(obj)
  if (length(features) > 1)
    df$order = colSums(obj@assays$RNA@counts[features,])
  else
    df$order = obj@assays$RNA@counts[features,]
  n = ncol(obj)
  
  i = 1
  df_others = data.frame()
  for (feature in features) {
    if (feature %in% rownames(obj)) {
      this_col = cols[i]
      # this_col = "blue"
      if (i == 1) {
        df$feature = obj@assays$RNA@data[feature,]
        m <- grDevices::colorRamp(c("lightgray", this_col))( (1:n)/n )
        df$col = colour_values(df$feature, palette = m)
      } else {
        this_values = obj@assays$RNA@data[feature,]
        m <- grDevices::colorRamp(c("lightgray", this_col))( (1:n)/n )
        this_cols = colour_values(this_values, palette = m)
        names(this_cols) = colnames(obj)
        
        feature_pos = colnames(obj)[which(obj@assays$RNA@counts[feature,] > 0 )]
        newRow = df[feature_pos,]
        newRow$feature = this_values[feature_pos]
        newRow$col = this_cols[feature_pos]
        df_others = rbind(df_others, newRow)
      }
      i = i + 1
    } else {
      print("Feature not found in rownames of object nor meta data.")
    }
  }
  
  df = rbind(df, df_others)
  
  df = df[order(df$order, decreasing = F),]
  print(head(df))
  p = ggplot(df, aes(UMAP_1, UMAP_2)) + geom_point(size = my.pt.size, colour = df$col, alpha = 0.4) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
  p = LabelClusters(p, "ident", unique(df$ident), repel = F, colour = "black")
  
  # Color title based on genes
  title_str = "<span>"
  for (i in 1:length(features)) {
    if (i == length(features)) {
      title_str = paste0(title_str, "<span style='color:", as.character(cols[i]), ";'>", features[i], "</span>")
    } else {
      title_str = paste0(title_str, "<span style='color:", as.character(cols[i]), ";'>", features[i], "</span> and ")
    }
    if (i %% 3 == 0) {
      title_str = paste0(title_str, "<br></br>")
    }
  }
  title_str = paste0(title_str, "</span>")
  p = p + labs(title = title_str) + theme( plot.title = element_markdown(lineheight = 1.1 ))
  
  return(p)
}

myFeaturePlot = function(obj, feature, cells.use = NULL, myslot = "data", alpha_vect = NULL, my.split.by = NULL, my.pt.size = 1, my.col.pal = NULL, my.title = NULL, na.blank = FALSE) {
  #' Make a FeaturePlot using ggplot. Having this gives me more control over the plot.
  #' 
  #' @param obj Seurat object
  #' @param feature feature to plot
  #' @param myslot slot to pull data from
  #' @param alpha_vect vector of values to use as an alpha parameter (values in same order as cells)
  #' @param my.split.by metadata to split plot by
  #' @param my.pt.size point size
  #' @param my.col.pal color palette
  #' @param my.title title to add (only if splitting plot)
  #' @param na.blank make cells that don't express the gene NA (that way they are colored gray)
  #' @return ggplot object
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  if (feature %in% rownames(obj)) {
    values = exp[feature,]
  } else if (feature %in% colnames(obj@meta.data)) {
    values = obj@meta.data[, feature]
  } else {
    print("Feature not found in rownames of object nor meta data.")
  }
  
  # Create a Dataframe with the UMAP coordinates, the expression values, and the Idents
  df = as.data.frame(obj@reductions$umap@cell.embeddings)
  df$value = values
  df$ident = Idents(obj)
  
  # Make the Cells that do not express the gene gray
  my.na.value = "grey90" # a constant used in the graphs
  if (na.blank) {
    if (feature %in% rownames(obj))
      df$value[which(obj@assays$RNA@counts[feature,] < 1)] = NA
    if (feature %in% colnames(obj@meta.data))
      df$value[which(df$value == 0)] = NA
  }
    
  
  onePlot = function(df, mymin, mymax) {
    if ( is.null(alpha_vect) ) {
      df = df[order(df$value, na.last = F, decreasing = F),]
      print(head(df))
      # p = ggplot(df, aes(UMAP_1, UMAP_2, col = value)) + geom_point(size = my.pt.size) + scale_color_gradient2(midpoint=(mymax-mymin)/2 + mymin, low = "blue", mid = "gold", high = "red", space = "Lab", limits=c(mymin, mymax)) + theme_classic() + theme(plot.title = element_text(hjust = 0.5))
      p = ggplot(df, aes(UMAP_1, UMAP_2, col = value)) + geom_point(size = my.pt.size) + scale_color_gradientn(limits = c(mymin, mymax), colors = c("lightgrey", "blue"), na.value=my.na.value) + theme_classic() + theme(plot.title = element_text(hjust = 0.5))
      if (!is.null(my.col.pal))
        p = p + scale_color_gradientn(colors = my.col.pal(50), na.value=my.na.value, limits = c(mymin, mymax))
    } else {
      rescale <- function(x, newMin, newMax) { (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin }
      df$total_counts = rescale(alpha_vect, 0.05, 0.6)[rownames(df)]
      df$total_counts = alpha_vect
      df = df[order(df$value),]
      print(head(df))
      p = ggplot(df, aes(UMAP_1, UMAP_2, col = value)) + geom_point(aes(size = total_counts, alpha = total_counts)) + scale_color_gradient2(midpoint=(mymax-mymin)/2 + mymin, low = "blue", mid = "gold", high = "red", space = "Lab", limits=c(mymin, mymax), na.value=my.na.value) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + scale_size(range=c(0.2, 2)) + scale_alpha(range=c(0.2, 1))
    }
    p = LabelClusters(p, "ident", unique(df$ident), repel = F, colour = "black")
    return(p)
  }
  
  # Find the min and max values (Used for scaling colors later on)
  if (! is.null(cells.use) )
    df = df[cells.use,]
  mymax = max(df$value[is.finite(df$value)], na.rm = T)
  mymin = min(df$value[is.finite(df$value)], na.rm = T)
  
  # Split Plot if Necessary
  if ( is.null(my.split.by) ) {
    p = onePlot(df, mymin, mymax)
    if ( is.null(my.title) )
      p = p + ggtitle(feature)
    else
      p = p + ggtitle(my.title)
  } else {
    p_list = list()
    split_unique = unique(obj@meta.data[c(my.split.by)][,c(my.split.by)])
    # split_unique = levels(bb@meta.data[c("sample")][,c("sample")])
    for (one_split in split_unique) {
      this_cells = rownames(df)[which( rownames(df) %in% colnames(obj)[which(obj@meta.data[c(my.split.by)] == one_split)] )]
      print(one_split)
      print(head(this_cells))
      p_list[[one_split]] = onePlot(df[this_cells,], mymin, mymax) + ggtitle(one_split)
    }
    p = plot_grid(plotlist=p_list, ncol = 2)
    
    # Add a Title if Given
    if (! is.null(my.title) ) {
      my.title.obj = ggdraw() + draw_label(my.title)
      p = plot_grid(my.title.obj, p, ncol = 1, rel_heights=c(0.1, 1))
    }
  } # more than one plot
  
  return(p)
}

gtfParser = function() {
  #' Junk function, just useful to have as a reference for iterating through a gtf.
  names = c()
  ids = c()
  for (i in 1:nrow(gtf)) {
    info = gtf[i,9]
    id_start = str_locate(pattern="gene_id", info)[2]
    id_stop = str_locate(pattern = ";", substr(info, id_start, 1000L))[1]
    name_start = str_locate(pattern = "gene_name", info)[2]
    name_stop = str_locate(pattern = ";", substr(info, name_start, 1000L))[1]
    
    id = substr(info, id_start+2, id_start+id_stop-2)
    name = substr(info, name_start+2, name_start+name_stop-2)
    names = c(names, name)
    ids = c(ids, id)
  }
  gtf$name = names
  gtf$id = ids
}

getDescription <- function(genes) {
  ensembl = useEnsembl(biomart="ensembl", dataset="mzebra_gene_ensembl")
  
  ensembl_description_all   <- unique(getBM(attributes=c('ensembl_gene_id', 'description'), filters = 'ensembl_gene_id', values = unique(genes), mart = ensembl))
  zfin_description_all      <- unique(getBM(attributes=c('zfin_id_symbol', 'description'), filters = 'zfin_id_symbol', values = unique(genes), mart = ensembl))
  hgnc_description_all      <- unique(getBM(attributes=c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = unique(genes), mart = ensembl))
  
  colnames(ensembl_description_all)   <- c('all_markers', 'all_markers_description')
  colnames(zfin_description_all)      <- c('all_markers', 'all_markers_description_2')
  colnames(hgnc_description_all)      <- c('all_markers', 'all_markers_description_3')
  
  ensembl_description_all = ensembl_description_all[!duplicated(ensembl_description_all$all_markers),]
  zfin_description_all    = zfin_description_all[!duplicated(zfin_description_all$all_markers),]
  hgnc_description_all    = hgnc_description_all[!duplicated(hgnc_description_all$all_markers),]
  
  total <- as.data.frame(genes)
  colnames(total) <- c("all_markers")
  try(total <- left_join(total, ensembl_description_all, by = c("all_markers")), silent = TRUE)
  try(total <- left_join(total, zfin_description_all, by = c("all_markers")), silent = TRUE)
  try(total <- left_join(total, hgnc_description_all, by = c("all_markers")), silent = TRUE)
  
  # Spacing
  try(total$all_markers_description[! is.na(total$all_markers_description)] <- total$all_markers_description[! is.na(total$all_markers_description)] + " ", silent = TRUE)
  try(total$all_markers_description_2[! is.na(total$all_markers_description_2)] <- total$all_markers_description_2[! is.na(total$all_markers_description_2)] + " ", silent = TRUE)
  try(total$all_markers_description_3[! is.na(total$all_markers_description_3)] <- total$all_markers_description_3[! is.na(total$all_markers_description_3)] + " ", silent = TRUE)
  
  try(total$all_markers_description[is.na(total$all_markers_description)] <- "", silent = TRUE)
  try(total$all_markers_description_2[is.na(total$all_markers_description_2)] <- "", silent = TRUE)
  try(total$all_markers_description_3[is.na(total$all_markers_description_3)] <- "", silent = TRUE)
  
  all_markers_description  = paste0(total$all_markers_description,  total$all_markers_description_2,  total$all_markers_description_3)
  # total$all_markers_description  <- all_markers_description
  # drops <- c("all_markers_description_2","cons_markers_description_2", "all_markers_description_3","cons_markers_description_3")
  # total <- total[ , !(names(total) %in% drops)]
  
  return(all_markers_description)
}

degWithHgncAndDescription = function(deg, gene_names) {
  #' Given a dataframe of degs, add a description and HGNC names for those genes
  #' 
  #' @param deg
  #' @return deg with 2 more columns. One for the description and one for the HGNC name.
  
  genes = deg$gene
  deg$description = getDescription(genes)
  gene_ind = which(colnames(deg) == "gene")
  deg$hgnc = deg$gene
  deg = deg[,c(1:gene_ind, ncol(deg), (gene_ind+1):(ncol(deg)-1))]
  deg = hgncMzebraInPlace(deg, gene_ind+1, gene_names, rm_na = F)
  return(deg)
}

cytoBINdeg = function(obj) {
  #' Find DEGs in each cytoBIN
  #' 
  #' @param obj Seurat object
  #' @return p1 barplot of Number of Cells in CytoBINs per cluster
  #' @return p2 barplot of Percent of Cells in CytoBINs per cluster
  #' @return deg DEGs per CytoBIN
  
  # Annotate each cell as eith Low, Medium, or High CytoBIN
  obj$cytoBIN = obj$cyto
  obj$cytoBIN[which(obj$cyto >= 0    & obj$cyto <= 0.33)] = "Low"
  obj$cytoBIN[which(obj$cyto > 0.33 & obj$cyto <= 0.66)] = "Medium"
  obj$cytoBIN[which(obj$cyto > 0.66 & obj$cyto <= 1   )] = "High"
  
  # See if the object identity is numeric
  clusters = sort(unique(as.vector(Idents(obj))))
  if (! any(is.na(as.numeric(Idents(obj))))) {
    clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  }
  
  # Find the number of cells in each cytoBIN for each cluster (or object identity)
  bin_df = data.frame()
  for (cluster in clusters) {
    this_cells = WhichCells(obj, idents = cluster)
    rows = data.frame(rep(cluster, 3), c("Low", "Medium", "High"),
                      c(length(this_cells[which(obj$cytoBIN[this_cells] == "Low")]),
                        length(this_cells[which(obj$cytoBIN[this_cells] == "Medium")]),
                        length(this_cells[which(obj$cytoBIN[this_cells] == "High")])))
    colnames(rows) = c("Cell_Type", "CytoBIN", "value")
    rows$pct = c(rows$value[1]/sum(rows$value), rows$value[2]/sum(rows$value), rows$value[3]/sum(rows$value))
    rows$pct = rows$pct*100
    bin_df = rbind(bin_df, rows)
  }
  
  # Plot the Number/Percent of Cells in Each CytoBIN for Each Cluster/Identity
  temp = rev(brewer.pal(11,"Spectral"))
  temp[6] = "gold" # this is what plotCytoTRACE uses
  bin_df$CytoBIN = factor(bin_df$CytoBIN, levels = c("High", "Medium", "Low"))
  bin_df$Cell_Type = factor(bin_df$Cell_Type, levels = clusters)
  p1 = ggplot(bin_df, aes(x=Cell_Type, y=value, fill=CytoBIN)) + geom_bar(stat="identity", position = "dodge") + scale_fill_manual(values=temp[c(10,6,2)]) + ggtitle("Number of Cells in CytoBINs per Cluster")
  p2 = ggplot(bin_df, aes(x=Cell_Type, y=pct, fill=CytoBIN)) + geom_bar(stat="identity", position = "dodge")   + scale_fill_manual(values=temp[c(10,6,2)]) + ggtitle("Percent of Cells in CytoBINs per Cluster")
  
  # Find DEGs per CytoBIN
  Idents(obj) = obj$cytoBIN
  deg = FindAllMarkers(obj, only.pos = T)
  deg = deg[which(deg$p_val_adj < 0.05),]
  
  return(list(p1, p2, deg))
}

cytoScoreByIdent = function(obj, my_idents = NULL, genes = NULL, by.gene = F, pt.alpha = 0.1) {
  #' Make a boxplot of cytoTRACE score by ident/cluster
  #' 
  #' @param obj Seurat object
  #' @param my_idents 
  #' @param genes plot the CytoTRACE score of cells positive for these genes (instead of plotting by identity)
  #' @return p boxplot of cytoTRACE score by ident/cluster
  
  # Color Palette
  temp = rev(brewer.pal(11,"Spectral"))
  temp[6] = "gold" # this is what plotCytoTRACE uses
  pal = colorRampPalette(temp)
  
  if (! is.null(genes) ) {
    # Plot by Gene Positive Cells
    cell_type_df = data.frame()
    for (gene in genes) {
      gene_pos = colnames(obj)[which(obj@assays$RNA@counts[gene,] > 0)]
      newRow = data.frame(CytoTRACE = as.numeric(as.vector(obj$cyto[gene_pos])), Cell_Type = gene)
      cell_type_df = rbind(cell_type_df, newRow)
    }
  } else {
    # Plot by Cell Types
    cell_type_df <- as.data.frame(cbind(as.numeric(as.vector(obj$cyto)), as.vector(Idents(obj))))
    colnames(cell_type_df) <- c("CytoTRACE", "Cell_Type")
    cell_type_df$CytoTRACE <- as.numeric(as.vector(cell_type_df$CytoTRACE))
    if (! is.null(my_idents) ) {
      # Plot by input cell types
      cell_type_df = cell_type_df[which(cell_type_df$Cell_Type %in% my_idents),]
      cell_type_df$Cell_Type <- factor(cell_type_df$Cell_Type, levels=my_idents)
    } else {
      # Plot by Cluster
      clusters = sort(unique(as.vector(Idents(obj))))
      if (! any(is.na(as.numeric(levels(Idents(obj)))))) {
        print("Can sort numerically")
        clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
      }
      cell_type_df$Cell_Type <- factor(cell_type_df$Cell_Type, levels=clusters)
    }
  }
  rownames(cell_type_df) <- NULL
  print(head(cell_type_df))
  print(length(which(is.na(as.numeric(as.vector(cell_type_df$CytoTRACE))))))
  cell_type_df$CytoTRACE = as.numeric(as.vector(cell_type_df$CytoTRACE))
  # p = ggplot(cell_type_df, aes(x=Cell_Type, y=CytoTRACE)) + geom_boxplot(alpha=0.8, aes(fill=as.character(as.numeric(as.vector(Cell_Type))%%2))) + ylim(0,1) + geom_jitter(position=position_dodge2(width=0.5), alpha=pt.alpha, aes(color = CytoTRACE)) + coord_flip() + theme_classic() + scale_color_gradientn(colors = pal(50)) + scale_fill_manual(values = c("white", "gray47")) + guides(color = F) + NoLegend()
  p = ggplot(cell_type_df, aes(x=Cell_Type, y=CytoTRACE)) + geom_boxplot(alpha=0.8, aes(fill=as.character(as.numeric(as.vector(Cell_Type))%%2))) + ylim(0,1) + geom_jitter(position=position_dodge2(width=0.5), alpha=pt.alpha, aes(color = CytoTRACE)) + coord_flip() + theme_classic() + scale_color_gradientn(colors = pal(50)) + scale_fill_manual(values = c("white", "gray47")) + guides(color = F) + NoLegend()
  if (! is.null(genes))
    p = p + ylab("Genes")
  return(p)
}

plotCyto = function(obj, cyto_feature = "cyto", png_name = "") {
  #' Plot CytoTRACE data
  #' 
  #' @param obj Seurat object to plot
  #' @param cyto_feature name of feature that holds CytoTRACE data (optional)
  #' @param png_name name of output png (optional)
  #' @return p UMAP of CytoTRACE score per cell
  temp = rev(brewer.pal(11,"Spectral"))
  temp[6] = "gold" # this is what plotCytoTRACE uses
  pal = colorRampPalette(temp)
  p = FeaturePlot(obj, features = "cyto", reduction = "umap", cols = pal(50), pt.size = 2, label=TRUE, order = TRUE) + guides(color=F)
  
  if (png_name != "") {
    png(width = 600, height = 600, file=png_name)
    print(p)
    dev.off()
  }
  
  return(p)
}

cellHeatmap = function(obj, markers, myslot="counts", non_zero=F) {
  #' Make a heamap of all the markers by all the cells
  #' 
  #' @param obj Seurat object
  #' @param markers vector gene markers
  #' @param myslot slot to pull expression data from
  #' @return p heatmap
  
  markers = markers[which(markers %in% rownames(obj))]
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  exp = exp[markers,]
  
  if (non_zero)
    exp = exp[,which(colSums(exp) > 0)]
  
  exp_melt = melt(as.matrix(exp))
  colnames(exp_melt) = c("Gene", "Cell", "value")
  
  p = ggplot(exp_melt, aes(Cell, Gene, fill=value)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, max(exp_melt$value))) + guides(color = FALSE) + theme_classic() + theme(axis.text.x = element_blank()) + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))
  return(p)
}

myDotPlot = function(obj, markers) {
  #' Dotplot of markers using Spectral colors. 
  #' And a bargraph of many genes are in each color for each cluster
  #' 
  #' @param obj Seurat Object
  #' @param markers vector of markers
  #' @return p DotPlot of markers
  # temp = rev(brewer.pal(11,"Spectral"))[c(rep(2,10), rep(10,10))]
  # temp = c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")
  temp = rev(c("#f94144", "#f3722c", "#f8961e", "#f9c74f", "#79bf43", "#43aa8b", "#577590"))
  # temp = rev(brewer.pal(11,"Spectral"))
  # temp[6] = "gold" # this is what plotCytoTRACE uses
  pal = colorRampPalette(temp)
  p = DotPlot(obj, features = markers, scale.by = "size", scale.min = 0) + scale_color_gradientn(colors = pal(50)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  # p = DotPlot(obj, features = markers, scale.by = "size", scale.min = 0) + scale_color_gradientn(colors = pal(3), values = c(1,0,-1)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  good_ind = which(complete.cases(p$data$avg.exp.scaled))
  min_avg = min(p$data$avg.exp.scaled[good_ind])
  max_avg = max(p$data$avg.exp.scaled[good_ind])
  third1 = ((max_avg-min_avg)*1/3) + min_avg
  third2 = ((max_avg-min_avg)*2/3) + min_avg
  
  p$data$third = p$data$avg.exp.scaled
  p$data$third[which(p$data$avg.exp.scaled >= min_avg & p$data$avg.exp.scaled < third1)]   = "Low"
  p$data$third[which(p$data$avg.exp.scaled >= third1  & p$data$avg.exp.scaled < third2)]   = "Medium"
  p$data$third[which(p$data$avg.exp.scaled >= third2  & p$data$avg.exp.scaled <= max_avg)] = "High"
  
  p$data$third = factor(p$data$third, levels = c("High", "Medium", "Low", "NaN"))
  # print(nrow(p$data[which(p$data$id == 0 & p$data$third == "High"),]))
  # print(p$data[which(p$data$id == 0 & p$data$third == "High"),])
  p1 = ggplot(p$data, aes(id, fill=third, color=third)) + geom_bar(position = position_dodge2(), alpha = 0.7) + scale_color_manual(values=c(temp[10], temp[6], temp[2], "gray")) + scale_fill_manual(values=c(temp[10], temp[6], temp[2], "gray")) + labs(color="Avg Expression", fill = "Avg Expression", x = "Cluster", y = "Number of Genes", title="Number of Genes in Each Avg Expression Bin Per Cluster") 
  # p2 = ggplot(p$data, aes(id, y=pct.exp, fill=third, color=third, group = third)) + geom_bar(stat="identity", alpha = 0.7) + scale_color_manual(values=c(temp[10], temp[6], temp[2], "gray")) + scale_fill_manual(values=c(temp[10], temp[6], temp[2], "gray")) + labs(color="Avg Expression", fill = "Avg Expression", x = "Cluster", y = "% Expressed", title="% Expressed for Each Gene in Each Avg Expression Bin Per Cluster") + facet_wrap(~ id)
  p2 = ggplot(p$data, aes(x=third, y=pct.exp, fill=third, color=third, group = third)) + geom_bar(stat="identity", alpha = 0.7) + scale_color_manual(values=c(temp[10], temp[6], temp[2], "gray")) + scale_fill_manual(values=c(temp[10], temp[6], temp[2], "gray")) + labs(color="Avg Expression", fill = "Avg Expression", x = "Cluster", y = "% Expressed", title="% Expressed for Each Gene in Each Avg Expression Bin Per Cluster") + facet_wrap(~ id)
  
  return(list(p, p1, p2))
}

markerLogFC = function(obj, markers, myslot="data") {
  #' Stacked barchart of the logFC of each marker for each cluster.
  #' Calculation for avg_logFC came from Seurat here:
  #' https://github.com/satijalab/seurat/issues/741
  #' 
  #' @param obj Seurat object
  #' @param markers vector of marker genes
  #' @param myslot slot of Seurat Object to pull expression data from
  #' @return list of 2
  #' 1) ggplot of a stacked barchart
  #' 2) dataframe used to build the barchart
  
  markers = markers[which(markers %in% rownames(obj))]
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  exp = exp[markers,]
  
  clusters = sort(unique(as.vector(Idents(obj))))
  if (! any(is.na(as.numeric(Idents(obj))))) {
    print("Can sort numerically")
    clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  }
  exp_df = data.frame()
  for (i in 1:length(clusters)) {
    cluster = clusters[i]
    cluster_cells <- WhichCells(obj, idents = cluster)
    other_cells = colnames(obj)[which(! colnames(obj) %in% cluster_cells)]
    
    if (length(markers) > 1) {
      cluster_mean_exp = rowMeans(expm1(exp[markers, cluster_cells]))
      other_mean_exp = rowMeans(expm1(exp[markers, other_cells]))
    } else {
      cluster_mean_exp = mean(expm1(exp[cluster_cells]))
      other_mean_exp = mean(expm1(exp[other_cells]))
    }
    
    avg_logFC = log(cluster_mean_exp + 1) - log(other_mean_exp + 1)
    exp_df = rbind(exp_df, cbind(rep(cluster, length(markers)), markers, avg_logFC))
  }
  print(head(exp_df))
  colnames(exp_df)[1] = "cluster"
  exp_df$avg_logFC = as.vector(as.numeric(exp_df$avg_logFC))
  exp_df$Regulation = exp_df$avg_logFC > 0
  exp_df$Regulation = plyr::revalue(as.character(exp_df$Regulation), replace = c("TRUE" = "Up", "FALSE" = "Down"))
  exp_df$cluster = factor(exp_df$cluster, levels=clusters)
  
  p = ggplot(exp_df, aes(x=cluster, y = avg_logFC, fill = Regulation, color = Regulation)) + geom_bar(stat="identity", alpha=0.6) + ylab("Average Log Fold Change") + geom_hline(yintercept=0) + ggtitle("Average LogFC for All Markers") + scale_fill_discrete(drop=FALSE, limits=levels(exp_df$Regulation)) + scale_color_discrete(drop=FALSE, limits=levels(exp_df$Regulation))
  
  return(list(p, exp_df))
}

markerHeatmap = function(obj, markers, myslot="scale.data") {
  #' Make a heatmap of all the markers by all the clusters
  #'
  #' @param obj Seurat object
  #' @param markers vector of marker genes
  #' @param myslot slot of Seurat Object to pull expression data from
  #' @return p heatmap
  
  markers = markers[which(markers %in% rownames(obj))]
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  exp = exp[markers,]
  
  # clusters = sort(unique(as.vector(Idents(obj))))
  clusters = levels(Idents(bb))
  if (! any(is.na(as.numeric(as.vector(Idents(obj)))))) {
    print("Can sort numerically")
    clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  }
  exp_df = data.frame()
  for (i in 1:length(clusters)) {
    cluster = clusters[i]
    cluster_cells <- WhichCells(obj, idents = cluster)
    cluster_mean_exp = rowMeans(exp[markers, cluster_cells])
    exp_df = rbind(exp_df, data.frame(rep(cluster, length(markers)), markers, cluster_mean_exp))
  }
  colnames(exp_df) = c("Cluster", "Gene", "Mean_Expression")
  
  # exp_df = as.data.table(exp_df)
  exp_df$Mean_Expression = as.numeric(as.vector(exp_df$Mean_Expression))
  exp_df$Gene = factor(exp_df$Gene, levels=markers)
  # exp_df[, fillScaled := scale(Mean_Expression), Gene] # scales value per gene
  
  # p = ggplot(exp_df, aes(Cluster, Gene)) + geom_tile(aes(fill=Mean_Expression)) + scale_fill_viridis(discrete=FALSE, breaks = c(min(exp_df$Mean_Expression), max(exp_df$Mean_Expression)), labels=c("Low", "High")) + guides(color = FALSE) + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + labs(fill = "Relative Expression") + scale_x_continuous(expand=c(0,0))
  p = ggplot(exp_df, aes(Gene, Cluster)) + geom_tile(aes(fill=Mean_Expression)) + scale_fill_viridis(discrete=FALSE) + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + labs(fill = "Scaled Expression")
  
  # Increase number of ticks
  x_labs = clusters
  x_labs[c(FALSE, TRUE)] = "" # R trickery -> selects every other element
  if (! any(is.na(as.numeric(as.vector(Idents(obj))))))
    p = p + scale_x_continuous(breaks=clusters, labels =x_labs, expand=c(0,0))
  
  return(p)
}

degDend = function(mat, genes, png_name, include_samples=c()) {
  #' Dendrogram of DEGs
  #' 
  #' @param mat matrix of counts data for all genes
  #' @param genes genes to use in the dendrogram
  #' @param include_samples only include samples in the dendrogram with these names
  
  if (length(include_samples) != 0) {
    mat = mat[,which(colnames(mat) %in% include_samples)]
  }
  mat = mat[which(rownames(mat) %in% genes),]
  
  # z = cpm(mat, normalized.lib.size=TRUE)
  # scaledata <- t(scale(t(z)))
  # scaledata <- scaledata[complete.cases(scaledata),]
  # mat = scaledata
  
  mat = log2(mat)
  mat <- mat[is.finite(rowSums(mat)),]
  
  Colors=rev(brewer.pal(11,"Spectral"))
  my_palette1 = colorRampPalette(Colors)(n = 299)
  my_palette <- colorRampPalette(c("#ff4b5c", "#FFFFFF", "#056674"))(n = 299)
  print("Plotting the dendrogram")
  png(png_name, width=500, height=500, res=100)
  heatmap.2(mat, scale = "none", dendrogram = "both", col = my_palette1, trace = "none", srtCol=45, symm=F,symkey=F,symbreaks=F, main="Log2 Expression")
  dev.off()
  p = "hi"
  
  return(p)
}

cellCycle = function(obj, isMzebra = T, isMouse=F, isNCBI = F, work = T) {
  #' Paint Seurat's Cell Cycle Annotation for each Cell
  #' 
  #' @param obj Seurat object
  #' @param isMzebra is the Seurat object from a cichlid?
  #' @param isMouse is the Seurat object from a mouse?
  #' @return cc_fact factor of cell cycle state for each cell (you can make this into metadata)
  
  if (work)
    rna_path = "~/research/"
  
  library("Seurat")
  s.genes   <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  if (isMzebra) {
    if (isMouse) {
      print("Mouse and Mzebra selected. Please choose only one organism (default is Mzebra).")
    }
    if (isNCBI) {
      gene_info = read.table(paste0(rna_path, "/all_research/gene_info.txt"), sep="\t", header = T, stringsAsFactors = F)
      s.genes   <- gene_info$mzebra[match(s.genes, gene_info$human)]
      g2m.genes <- gene_info$mzebra[match(g2m.genes, gene_info$human)]
      s.genes = s.genes[which(! is.na(s.genes) )]
      g2m.genes = s.genes[which(! is.na(g2m.genes) )]
    } else {
      s.genes   <- convertHgncDataFrameToMzebra(data.frame(s.genes), gene_column = 1, gene_names = rownames(obj), na.rm = T, return_vect = T)
      g2m.genes <- convertHgncDataFrameToMzebra(data.frame(g2m.genes), gene_column = 1, gene_names = rownames(obj), na.rm = T, return_vect = T)
    }
  } else if (isMouse) {
    s.genes = str_to_title(s.genes)
    g2m.genes = str_to_title(g2m.genes)
  }
  
  obj = CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  # p = DimPlot(obj, label = T, order = T)
  
  return(obj$Phase)
}

numCellExpMarker = function(obj, markers, myslot="data", thresh=0) {
  #' Find the Number of Cells Expressing the Markers in Comparison to A List of Random
  #' Genes of the Same Length At a Threshold.
  #' 
  #' @param obj Seurat Object
  #' @param markers vector of marker genes
  #' @param myslot slot to pull data from (either "counts" or "data")
  #' @return p histogram of cells expressing the markers vs random list
  
  # Generate a list of random genes of equal size as the markers
  ran = sample(rownames(obj), length(markers), replace = F)
  
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  clusters = sort(unique(as.numeric(as.vector(obj$seurat_clusters))))
  num_cells = data.frame()
  for (gene in ran) {
    num_cells <- rbind(num_cells, t(c("Random", length(which(as.vector(exp[gene,]) != 0)))))
  }
  for (gene in markers) {
    num_cells <- rbind(num_cells, t(c("Real", length(which(as.vector(exp[gene,]) != 0)))))
  }
  colnames(num_cells) <- c("type", "num_cells")
  num_cells$num_cells <- as.numeric(as.vector(num_cells$num_cells)) + 0.1
  p = ggplot(num_cells, aes(x=num_cells, fill = type, color=type)) + geom_histogram(stat = "bin", alpha=0.5, position="identity") + ggtitle("# of Cells Expressing a Random List of Genes vs List of Markers") + xlab("Number of Cells") + ylab("Number of Genes")
  
  return(p)
}

markerExpPerCell = function(obj, markers, myslot="data", n_markers=F) {
  #' Paint the average Expression of Markers Per Cell (UMAP)
  #' 
  #' @param obj Seurat Object
  #' @param markers vector of markers
  #' @param myslot slot to pull data from (either "counts" or "data")
  #' @param n_markers plot the number of markers expressed per cell instead?
  #' @return p ggplot of expression in every cell
  
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  exp = exp[markers,]
  if (n_markers) {
    exp[which(exp > 0)] = 1
    obj$value = colSums(exp)
    title_str = "Number of Markers Expressed per Cell"
  } else {
    obj$value = colSums(exp)
    title_str = "Expression of Markers per Cell"
  }
  
  p = FeaturePlot(obj, features = "value", order = TRUE, label = TRUE, pt.size = 1.5) + labs(title = title_str)
  return(p)
}

markerCellPerCluster = function(obj, markers, correct=T) {
  #' For each marker find the number of positive cells in that cluster.
  #' Plot the cummulative number and by default correct for cluster size and background expression.
  #' 
  #' @param obj Seurat object
  #' @param markers vector of markers
  #' @param correct correct for cluster size and background expression?
  #' @return p barplot of cells
  
  exp = GetAssayData(obj, assay = "RNA", slot="counts")
  exp[which(exp > 0)] = 1
  
  title_str = if_else(correct, "Normalized", "")
  title_str = paste(title_str, "Total Number of Cells per Cluster Expressing a Marker")
  
  per_cluster_df = data.frame()
  clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  for (cluster in clusters) {
    print(cluster)
    cluster_cells = WhichCells(obj, idents = cluster)
    if (length(markers) > 1)
      num_cells = sum(colSums(exp[markers, cluster_cells]))
    else
      num_cells = sum(exp[markers, cluster_cells])
    if (correct)
      num_cells = num_cells / sum(colSums(exp[, cluster_cells]))
    print(num_cells)
    per_cluster_df = rbind(per_cluster_df, data.frame(cluster, num_cells))
  }
  colnames(per_cluster_df) = c("cluster", "num_cells")
  per_cluster_df$cluster = factor(per_cluster_df$cluster, levels = clusters)
  p = ggplot(per_cluster_df, aes(cluster, num_cells)) + geom_bar(stat = "identity") + ggtitle(title_str)
  
  return(p)
}

markerExpPerCellPerClusterBVC = function(obj, markers, myslot="counts", n_markers=T, correct=F) {
  #' *This is the new function that tests behave vs control.*
  #' Boxplot of Marker Expression Per Cell Per Cluster
  #' 
  #' @param obj Seurat Object
  #' @param markers vector of markers
  #' @param n_markers count the number of markers per cell?
  #' @param correct correct for background expresssion?
  #' @param myslot slot to pull data from (either "counts" or "data")
  #' @return list of 3:
  #' 1) boxplot with jitter points of marker expression per cell per cluster
  #' 2) effect size
  #' 3) dataframe with columns: 
  #'    cluster, magnitude of Cohen's d, Cohen's d, Upper bound on Cohen's confidence interval, Lower bound on Cohen's confidence interval, p-value from ztest, bonferroni corrected p-value
  
  title_str = "Avreage"
  title_str = paste(title_str, if_else(myslot=="counts", "", "Normalized"))
  title_str = paste(title_str, "Expression of Marker Genes per Cell per Cluster")
  
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  if (n_markers) {
    title_str = "Number of Markers Genes Expressed per Cell per Cluster"
    exp[which(exp > 0)] = 1
  }
  title_str = paste(title_str, if_else(correct, "- Corrected", ""))
  
  # Idents(obj) = obj$seurat_clusters
  per_cluster_df = data.frame()
  d_df = data.frame()
  clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  conds = unique(obj$cond)
  for (cluster in clusters) {
    cluster_cells <- WhichCells(obj, idents = cluster)
    for (cond in conds) {
      cond_cells = colnames(obj)[which(obj$cond == cond)]
      cond_cluster_cells = cluster_cells[which(cluster_cells %in% cond_cells)]
      avg_cluster_exp = colSums(exp[markers, cond_cluster_cells])
      if (correct)
        avg_cluster_exp = avg_cluster_exp/colSums(exp[, cond_cluster_cells])
      
      # EDIT
      other_cells = colnames(obj)[which(! colnames(obj) %in% cluster_cells)]
      other_avg = colSums(exp[markers, other_cells])
      if (correct)
        other_avg = other_avg/colSums(exp[,other_cells])
      
      # Z-test
      p = z.test(avg_cluster_exp, other_avg, sigma.x = sd(avg_cluster_exp), sigma.y = sd(other_avg), alternative = "greater")$p.value
      
      # Cohen's d
      all_exp = c(avg_cluster_exp, other_avg)
      test= effsize::cohen.d(all_exp, c(rep("cluster", length(avg_cluster_exp)), 
                                        rep("other",   length(other_avg))))
      
      d=test$estimate
      up=test$conf.int[2]
      down = test$conf.int[1]
      mag=test$magnitude
      mag_pos=mag
      mag_pos[which(d < 0)] = "negligible"
      # WILCOX
      # p_val = wilcox.test(avg_cluster_exp, other_avg, alternative = "greater")$p.value
      
      d_df = rbind(d_df, data.frame(cluster, cond, mag, mag_pos, d, up, down, p))
      per_cluster_df = rbind(per_cluster_df, data.frame(cluster, cond, avg_cluster_exp, p, mag_pos))
    } # end condition for
    
  } # end clusters for
  
  # mycols = c("#01c5c4", "#b8de6f", "#f1e189", "#f39233")
  my_levels = apply(rev(expand.grid(unique(obj$cond), clusters)), 1, paste, collapse=" ")
  my_levels = trimws(my_levels)
  xtext_cols = rep(c("red", "blue"), length(clusters))
  d_df$cluster_cond = paste(d_df$cluster, d_df$cond)
  d_df$p_val_adj = p.adjust(d_df$p, method = "bonferroni")
  d_df$cluster_cond = factor(d_df$cluster_cond, levels = my_levels)
  p1  = ggplot(d_df, aes(cluster_cond, d, color=mag, fill = mag)) + geom_pointrange(aes(ymin = down, ymax = up)) + ylab("Cohen's d") + ggtitle("Effect Size in Clusters") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols))
  
  colnames(per_cluster_df) <- c("cluster", "cond", "avg_cluster_exp", "p", "mag_pos")
  per_cluster_df$avg_cluster_exp = as.numeric(as.vector(per_cluster_df$avg_cluster_exp))
  per_cluster_df$cluster_cond = paste(per_cluster_df$cluster, per_cluster_df$cond)
  per_cluster_df$p_val_adj = d_df$p_val_adj[match(per_cluster_df$cluster_cond, d_df$cluster_cond)]
  per_cluster_df$star = ifelse(per_cluster_df$p_val_adj < 0.001, "***", 
                               ifelse(per_cluster_df$p_val_adj < 0.01, "**", 
                                      ifelse(per_cluster_df$p_val_adj < 0.05, "*", "")))
  per_cluster_df$cluster = factor(per_cluster_df$cluster, levels = clusters)
  per_cluster_df$cluster_cond = factor(per_cluster_df$cluster_cond, levels = my_levels)
  per_cluster_df$mag_pos = factor(per_cluster_df$mag_pos, levels = c("negligible", "small", "medium", "large"))
  p = ggplot(per_cluster_df, aes(cluster_cond, avg_cluster_exp, fill=mag_pos, color=mag_pos)) + geom_text(aes(x= cluster_cond, y = Inf, vjust = 1, label = star)) + geom_boxplot(alpha=0.6) +  geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.05) + xlab("Cluster") + ylab("") + ggtitle(title_str) + scale_color_viridis(discrete = T,drop=TRUE, limits = levels(per_cluster_df$mag_pos), name = "Effect Size") + scale_fill_viridis(discrete = T, drop=TRUE, limits = levels(per_cluster_df$mag_pos), name = "Effect Size") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols))
  
  print(head(per_cluster_df[which(per_cluster_df$cluster == 0),]))
  
  return(list(p, p1, d_df))
}

separateEnrichTest = function(obj, markers) {
  #' Experimental Enrichment Test. Plot each gene separately
  exp = GetAssayData(obj, assay = "RNA", slot="counts")
  exp = exp[markers,]
  exp[which(exp > 0)] = 1
  
  per_cluster_df = data.frame()
  per_gene_per_cluster_df = data.frame()
  clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  for (cluster in clusters) {
    cluster_cells <- WhichCells(obj, idents = cluster)
    other_cells = colnames(obj)[which(! colnames(obj) %in% cluster_cells )]
    cluster_ps = c()
    for (gene in markers) {
      this_exp = exp[gene,]
      this_exp_correct = exp[gene,] / obj$nFeature_RNA
      c_exp_cor = this_exp_correct[cluster_cells]  # cluster expression corrected
      o_exp_cor = this_exp_correct[other_cells]  # other expression corrected
      p = z.test(c_exp_cor, o_exp_cor, sigma.x = sd(c_exp_cor), sigma.y = sd(o_exp_cor), alternative = "greater")$p.value
      cluster_ps = c(cluster_ps, p)
      per_gene_per_cluster_df = rbind(per_gene_per_cluster_df, t(c(cluster, gene, p)))
    }
    composite_p = sumz(cluster_ps)$p[1,1]
    per_cluster_df = rbind(per_cluster_df, t(c(cluster, cluster_ps, composite_p)))
  }
  colnames(per_cluster_df) = c("cluster", markers, "composite_p")
  colnames(per_gene_per_cluster_df) = c("cluster", "gene", "p")
  per_cluster_df$p_adj = p.adjust(per_cluster_df$composite_p, method = "bonferroni")
  
  return(list(per_gene_per_cluster_df, per_cluster_df))
}

markerExpPerCellPerClusterQuickGeneX = function(obj, markers, geneX) {
  #' Same as markerExpPerCellPerClusterQuick but on GeneX+ vs GeneX-
  # Defaults from original function
  n_markers = T
  correct = T
  myslot = "counts"
  
  title_str = "Average"
  title_str = paste(title_str, if_else(myslot=="counts", "", "Normalized"))
  title_str = paste(title_str, "Expression of Marker Genes per Cell per Cluster")
  
  obj$geneX = "Neg"
  obj$geneX[which(obj@assays$RNA@counts[geneX,] > 0)] = "Pos"
  if (length(which(obj$geneX == "Pos")) < 3)
    return(NA)
  Idents(obj) = obj$geneX
  
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  exp = exp[markers,]
  exp[which(exp > 0)] = 1
  d_df = data.frame()
  clusters = c("Pos", "Neg")
  
  for (cluster in clusters) {
    cluster_cells <- WhichCells(obj, idents = cluster)
    num_cells = length(cluster_cells)
    
    avg_cluster_exp = colSums(exp[markers, cluster_cells])
    avg_cluster_exp = avg_cluster_exp/obj$nFeature_RNA[cluster_cells]
    other_cells = colnames(obj)[which(! colnames(obj) %in% cluster_cells)]
    other_avg = colSums(exp[markers, other_cells])
    other_avg = other_avg/obj$nFeature_RNA[other_cells]
    
    p = z.test(avg_cluster_exp, other_avg, sigma.x = sd(avg_cluster_exp), sigma.y = sd(other_avg), alternative = "greater")$p.value
    
    # Cohen's d
    all_exp = c(avg_cluster_exp, other_avg)
    test= effsize::cohen.d(all_exp, c(rep("cluster", length(avg_cluster_exp)),
                                      rep("other",   length(other_avg))))
    
    d=test$estimate
    up=test$conf.int[2]
    down = test$conf.int[1]
    mag=test$magnitude
    mag_pos=mag
    mag_pos[which(d < 0)] = "negligible"
    
    d_df = rbind(d_df, data.frame(cluster, mag, mag_pos, d, up, down, num_cells, p))
  }
  d_df$cluster = factor(d_df$cluster, levels = clusters)
  
  return(d_df)
}

markerExpPerCellPerClusterQuick = function(obj, markers) {
  #' Same as markerExpPerCellPerCluster using the default settings and using some shortcuts
  # Defaults from original function
  n_markers = T
  correct = T
  myslot = "counts"
  
  title_str = "Average"
  title_str = paste(title_str, if_else(myslot=="counts", "", "Normalized"))
  title_str = paste(title_str, "Expression of Marker Genes per Cell per Cluster")

  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  exp = exp[markers,]
  exp[which(exp > 0)] = 1
  per_cluster_df = data.frame()
  d_df = data.frame()
  clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  
  for (cluster in clusters) {
    cluster_cells <- WhichCells(obj, idents = cluster)

    avg_cluster_exp = colSums(exp[markers, cluster_cells])
    avg_cluster_exp = avg_cluster_exp/obj$nFeature_RNA[cluster_cells]
    other_cells = colnames(obj)[which(! colnames(obj) %in% cluster_cells)]
    other_avg = colSums(exp[markers, other_cells])
    other_avg = other_avg/obj$nFeature_RNA[other_cells]

    p = z.test(avg_cluster_exp, other_avg, sigma.x = sd(avg_cluster_exp), sigma.y = sd(other_avg), alternative = "greater")$p.value

    # Cohen's d
    all_exp = c(avg_cluster_exp, other_avg)
    test= effsize::cohen.d(all_exp, c(rep("cluster", length(avg_cluster_exp)),
                                      rep("other",   length(other_avg))))

    d=test$estimate
    up=test$conf.int[2]
    down = test$conf.int[1]
    mag=test$magnitude
    mag_pos=mag
    mag_pos[which(d < 0)] = "negligible"

    d_df = rbind(d_df, data.frame(cluster, mag, mag_pos, d, up, down, p))
    per_cluster_df = rbind(per_cluster_df, data.frame(cluster, avg_cluster_exp, p, mag_pos))
  }
  
  d_df$p_val_adj = p.adjust(d_df$p, method = "bonferroni")
  d_df$cluster = factor(d_df$cluster, levels = clusters)
  p1  = ggplot(d_df, aes(cluster, d, color=mag, fill = mag)) + geom_pointrange(aes(ymin = down, ymax = up)) + ylab("Cohen's d") + ggtitle("Effect Size in Clusters")

  colnames(per_cluster_df) <- c("cluster", "avg_cluster_exp", "p", "mag_pos")
  per_cluster_df$avg_cluster_exp = as.numeric(as.vector(per_cluster_df$avg_cluster_exp))
  per_cluster_df$p_val_adj = unlist(sapply(1:length(clusters), function(x) rep(d_df$p_val_adj[which(d_df$cluster == clusters[x])], length(which(per_cluster_df$cluster == clusters[x]))) ))
  per_cluster_df$star = ifelse(per_cluster_df$p_val_adj < 0.001, "***",
                               ifelse(per_cluster_df$p_val_adj < 0.01, "**",
                                      ifelse(per_cluster_df$p_val_adj < 0.05, "*", "")))
  per_cluster_df$cluster = factor(per_cluster_df$cluster, levels = clusters)
  per_cluster_df$mag_pos = factor(per_cluster_df$mag_pos, levels = c("negligible", "small", "medium", "large"))
  p = ggplot(per_cluster_df, aes(cluster, avg_cluster_exp, fill=mag_pos, color=mag_pos)) + geom_text(aes(x= cluster, y = Inf, vjust = 1, label = star)) + geom_boxplot(alpha=0.6) +  geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.05) + xlab("Cluster") + ylab("") + ggtitle(title_str) + scale_color_viridis(discrete = T,drop=TRUE, limits = levels(per_cluster_df$mag_pos), name = "Effect Size") + scale_fill_viridis(discrete = T, drop=TRUE, limits = levels(per_cluster_df$mag_pos), name = "Effect Size")

  print(head(per_cluster_df[which(per_cluster_df$cluster == 0),]))

  return(list(p, p1, d_df))
}

markerExpPerCellPerCluster = function(obj, markers, myslot="counts", n_markers=T, correct=T) {
  #' Boxplot of Marker Expression Per Cell Per Cluster
  #' 
  #' @param obj Seurat Object
  #' @param markers vector of markers
  #' @param n_markers count the number of markers per cell?
  #' @param correct correct for background expresssion?
  #' @param myslot slot to pull data from (either "counts" or "data")
  #' @return list of 3:
  #' 1) boxplot with jitter points of marker expression per cell per cluster
  #' 2) effect size
  #' 3) dataframe with columns: 
  #'    cluster, magnitude of Cohen's d, Cohen's d, Upper bound on Cohen's confidence interval, Lower bound on Cohen's confidence interval, p-value from ztest, bonferroni corrected p-value
  
  title_str = "Average"
  title_str = paste(title_str, if_else(myslot=="counts", "", "Normalized"))
  title_str = paste(title_str, "Expression of Marker Genes per Cell per Cluster")
  
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  if (n_markers) {
    title_str = "Number of Markers Genes Expressed per Cell per Cluster"
    exp[which(exp > 0)] = 1
  }
  title_str = paste(title_str, if_else(correct, "- Corrected", ""))
  
  # Idents(obj) = obj$seurat_clusters
  per_cluster_df = data.frame()
  d_df = data.frame()
  clusters = sort(unique(as.numeric(as.vector(Idents(obj)))))
  for (cluster in clusters) {
    cluster_cells <- WhichCells(obj, idents = cluster)
    
    # all_cluster_exp = colSums(exp[markers,]) / obj$nCount_RNA
    # all_cluster_exp_scale = scale(all_cluster_exp)
    # names(all_cluster_exp_scale) = colnames(obj)
    # avg_cluster_exp = all_cluster_exp_scale[which(names(all_cluster_exp_scale) %in% cluster_cells)]
    # other_avg = all_cluster_exp_scale[which(! names(all_cluster_exp_scale) %in% cluster_cells)]
    
    
    avg_cluster_exp = colSums(exp[markers, cluster_cells])
    if (correct)
      avg_cluster_exp = avg_cluster_exp/colSums(exp[,cluster_cells])

    # EDIT
    other_cells = colnames(obj)[which(! colnames(obj) %in% cluster_cells)]
    other_avg = colSums(exp[markers, other_cells])
    if (correct)
      other_avg = other_avg/colSums(exp[,other_cells])
    
    # Z-test
    p = z.test(avg_cluster_exp, other_avg, sigma.x = sd(avg_cluster_exp), sigma.y = sd(other_avg), alternative = "greater")$p.value
    
    # Cohen's d
    all_exp = c(avg_cluster_exp, other_avg)
    test= effsize::cohen.d(all_exp, c(rep("cluster", length(avg_cluster_exp)), 
                                      rep("other",   length(other_avg))))
    
    d=test$estimate
    up=test$conf.int[2]
    down = test$conf.int[1]
    mag=test$magnitude
    mag_pos=mag
    mag_pos[which(d < 0)] = "negligible"
    # WILCOX
    # p_val = wilcox.test(avg_cluster_exp, other_avg, alternative = "greater")$p.value
    
    d_df = rbind(d_df, data.frame(cluster, mag, mag_pos, d, up, down, p))
    per_cluster_df = rbind(per_cluster_df, data.frame(cluster, avg_cluster_exp, p, mag_pos))
  }
  
  # mycols = c("#01c5c4", "#b8de6f", "#f1e189", "#f39233")
  
  d_df$p_val_adj = p.adjust(d_df$p, method = "bonferroni")
  d_df$cluster = factor(d_df$cluster, levels = clusters)
  p1  = ggplot(d_df, aes(cluster, d, color=mag, fill = mag)) + geom_pointrange(aes(ymin = down, ymax = up)) + ylab("Cohen's d") + ggtitle("Effect Size in Clusters")
  
  colnames(per_cluster_df) <- c("cluster", "avg_cluster_exp", "p", "mag_pos")
  per_cluster_df$avg_cluster_exp = as.numeric(as.vector(per_cluster_df$avg_cluster_exp))
  per_cluster_df$p_val_adj = unlist(sapply(1:length(clusters), function(x) rep(d_df$p_val_adj[which(d_df$cluster == clusters[x])], length(which(per_cluster_df$cluster == clusters[x]))) ))
  per_cluster_df$star = ifelse(per_cluster_df$p_val_adj < 0.001, "***", 
                               ifelse(per_cluster_df$p_val_adj < 0.01, "**", 
                                      ifelse(per_cluster_df$p_val_adj < 0.05, "*", "")))
  per_cluster_df$cluster = factor(per_cluster_df$cluster, levels = clusters)
  per_cluster_df$mag_pos = factor(per_cluster_df$mag_pos, levels = c("negligible", "small", "medium", "large"))
  p = ggplot(per_cluster_df, aes(cluster, avg_cluster_exp, fill=mag_pos, color=mag_pos)) + geom_text(aes(x= cluster, y = Inf, vjust = 1, label = star)) + geom_boxplot(alpha=0.6) +  geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.05) + xlab("Cluster") + ylab("") + ggtitle(title_str) + scale_color_viridis(discrete = T,drop=TRUE, limits = levels(per_cluster_df$mag_pos), name = "Effect Size") + scale_fill_viridis(discrete = T, drop=TRUE, limits = levels(per_cluster_df$mag_pos), name = "Effect Size")
  
  print(head(per_cluster_df[which(per_cluster_df$cluster == 0),]))
  # ANOVA
  # per_cluster_df$cluster[which(per_cluster_df$cluster != 0)] = "all"
  # test = aov(avg_cluster_exp ~ cluster, data = per_cluster_df)
  # print(summary(test))
  # print(TukeyHSD(test))
  
  return(list(p, p1, d_df))
}

numMarkerTransPerCluster = function(obj, markers, myslot="counts", correct = T) {
  #' Find the Number of Marker Transcripts Per Cluster.
  #' 
  #' @param obj Seurat Object
  #' @param markers vector of markers
  #' @param myslot slot to pull data from (either "counts" or "data")
  #' @return p barchart of expression of marker per cluster. Each marker is a color.
  
  title_str = "Number of"
  title_str = paste(title_str, if_else(myslot=="counts", "Transcripts", "Normalized Transcripts"))
  title_str = paste(title_str, "of Marker Genes")
  
  exp = GetAssayData(obj, assay = "RNA", slot=myslot)
  # exp = exp[markers,]
  
  Idents(obj) = obj$seurat_clusters
  per_cluster_df = data.frame()
  correct_df = data.frame()
  clusters = sort(unique(as.numeric(as.vector(obj$seurat_clusters))))
  for (cluster in clusters) {
    cluster_cells <- WhichCells(obj, idents = cluster)
    cluster_sums = rowSums(exp[markers, cluster_cells])
    per_cluster_df = rbind(per_cluster_df, data.frame(markers, rep(cluster, length(markers)), cluster_sums))
    
    if (correct) {
      other_genes = rownames(obj)[which(! rownames(obj) %in% valid_genes)]
      cluster_sums = sum(cluster_sums)/sum(rowSums(exp[other_genes, cluster_cells]))
      correct_df = rbind(correct_df, data.frame(cluster, cluster_sums))
    }
  } # end cluster for
  
  colnames(per_cluster_df) <- c("gene", "cluster", "cluster_sum")
  p = ggplot(per_cluster_df, aes(cluster, as.numeric(as.vector(cluster_sum)), fill=gene)) + geom_bar(stat = "identity") + NoLegend() + xlab("Cluster") + ylab("") + ggtitle(title_str)
  
  if (correct) {
    colnames(correct_df) = c("cluster", "cluster_sum")
    my_mean = mean(correct_df$cluster_sum)
    print(my_mean)
    correct_df$cluster = factor(correct_df$cluster, levels = levels(Idents(obj)))
    p = ggplot(correct_df, aes(cluster, cluster_sum)) + geom_bar(stat = "identity") + xlab("Cluster") + ylab("Proportion of Transcripts Belonging to Markers") + geom_text(aes(y=my_mean, label=paste("Mean =", my_mean)), x = Inf, hjust=1, vjust=0, color = "red") + geom_hline(aes(yintercept = my_mean), color = "red")
  }
  
  return(p)
}

# This function produces the warning about column: 'cluster'
markersInDEGUpDown = function(deg, markers, gene_column = "gene", pct = F) {
  #' Search for markers that are in the up and downregulated DEGs
  #' 
  #' @param deg dataframe of DEGs (default output of FindAllMarkers)
  #' @param markers vector of markers
  #' @param gene_column name of the column that contains the genes
  #' @param pct plot pct of DEGs per cluster?
  #' @return results dataframe w/ col for marker and col for cluster the marker was found in
  #' @return p barchart of number of markers in cluster's degs
  
  clusters = unique(deg$cluster)
  if (! any(is.na(as.numeric(as.vector(clusters))))) {
    clusters = sort(as.numeric(as.vector(clusters)))
  }
  
  deg$isPos = deg$avg_logFC > 0
  
  gene_column_n = which(colnames(deg) == gene_column) # this is number
  results <- data.frame()
  for (gene in markers) {
    this_rows = deg[which(deg[,gene_column_n] == gene),]
    results = rbind(results, data.frame(rep(gene, nrow(this_rows)), this_rows$cluster, this_rows$avg_logFC, this_rows$isPos))
  } # end for
  if (nrow(results) > 0) {
    colnames(results) = c("gene", "cluster", "avg_logFC", "isPos")
    results$cluster = factor(results$cluster, levels = clusters)
    
    pos_res = table(results[which(results$isPos == T),c("gene", "cluster")])
    neg_res = table(results[which(results$isPos == F),c("gene", "cluster")])
    df = data.frame(up = colSums(pos_res), down = -colSums(neg_res), cluster = colnames(pos_res))
    df$cluster = factor(df$cluster, levels=clusters)
    colnames(df) = c("Up", "Down", "cluster")
    df2 = df %>% pivot_longer(c("Up", "Down"), names_to = "Regulation", values_to = "value")
    df2 = as.data.frame(df2)
    
    p = ggplot(df2, aes(x=cluster, y=value, fill=Regulation)) + geom_bar(stat="identity", alpha=0.8) + scale_x_discrete(drop=FALSE) + xlab("Cluster") + ylab("Number of Markers in Cluster DEGs") + geom_hline(yintercept=0)
    
    df3 = data.frame(cluster = clusters)
    df3$avg_logFC = sapply(1:length(clusters), function(x) sum(results$avg_logFC[which(results$cluster == clusters[[x]])]))
    df3$isPos = df3$avg_logFC > 0
    # p1 = ggplot(df3, aes(x=cluster, y = avg_logFC, fill = isPos)) + geom_bar(stat="identity") + guides(fill = F) + ylab("Total Sum of avg_logFC of Markers in DEGs") + xlab("Cluster") + geom_hline(yintercept=0)
    # p1 = ggplot(results, aes(x=cluster, y = avg_logFC, fill = isPos, color = isPos)) + geom_boxplot(alpha=0.6) + geom_jitter(position=position_dodge2(width=0.6), alpha=0.3)
    results$isPos = plyr::revalue(as.character(results$isPos), replace = c("TRUE" = "Up", "FALSE" = "Down"))
    colnames(results)[ncol(results)] = "Regulation"
    results$Regulation = factor(results$Regulation, levels = c("Down", "Up"))
    p1 = ggplot(results, aes(x=cluster, y = avg_logFC, fill = Regulation, color = Regulation)) + geom_bar(stat="identity", alpha=0.6) + ylab("Average Log Fold Change") + geom_hline(yintercept=0) + ggtitle("Average LogFC for Markers in DEGs") + scale_fill_discrete(drop=FALSE, limits=levels(results$Regulation)) + scale_color_discrete(drop=FALSE, limits=levels(results$Regulation))
    
    if (pct) {
      df2$pct = sapply(1:nrow(df2), function(x) df2$value[x]/nrow(deg[which(deg$cluster == df2$cluster[x]),]) * 100)
      p = ggplot(df2, aes(x=cluster, y=pct, fill=Regulation)) + geom_bar(stat="identity", alpha=0.8) + scale_x_discrete(drop=FALSE) + xlab("Cluster") + ylab("% of Cluster DEGs in Markers") + geom_hline(yintercept = 0)
    }
  } else {
    p = NULL
    p1 = NULL
    print("No Markers Found in the DEGs.")
  } # end else
  
  return(list(results, p, p1))
}

findMarkersInDEG = function(deg, markers, gene_column = "gene", pct = F) {
  #' Search for markers that are in the only upregulated DEGs
  #' 
  #' @param deg dataframe of DEGs (default output of FindAllMarkers)
  #' @param markers vector of markers
  #' @param gene_column name of the column that contains the genes
  #' @param pct plot pct of DEGs per cluster?
  #' @return results dataframe w/ col for marker and col for cluster the marker was found in
  #' @return p barchart of number of markers in cluster's degs
  
  clusters = unique(deg$cluster)
  if (! any(is.na(as.numeric(as.vector(clusters))))) {
    clusters = sort(as.numeric(as.vector(clusters)))
  }
  
  gene_column_n = which(colnames(deg) == gene_column) # this is number
  results <- data.frame()
  for (gene in markers) {
    this_rows = deg[which(deg[,gene_column_n] == gene),]
    results = rbind(results, data.frame(rep(gene, nrow(this_rows)), this_rows$cluster))
  } # end for
  if (nrow(results) > 0) {
    colnames(results) = c("gene", "cluster")
    results$cluster = factor(results$cluster, levels = clusters)
    p = ggplot(results, aes(cluster)) + geom_bar() + scale_x_discrete(drop=FALSE) + xlab("Cluster") + ylab("Number of Markers in Cluster DEGs")
    
    if (pct) {
      res_table = table(results)
      pct_df = data.frame(counts = colSums(res_table), cluster = colnames(res_table))
      pct_df$pct = sapply(1:length(clusters), function(x) pct_df$counts[which(pct_df$cluster == clusters[x])]/nrow(deg[which(deg$cluster == clusters[x]),]) * 100)
      p = ggplot(pct_df, aes(x=cluster, y=pct)) + geom_bar(stat = "identity") + scale_x_discrete(drop=FALSE) + xlab("Cluster") + ylab("% of Cluster DEGs in Markers")
    }
  } else {
    p = NULL
    print("No Markers Found in the DEGs.")
  } # end else
  
  return(list(results, p))
}

paintMarkers = function(obj, markers, filepath, mywidth=1200, myheight=500, myptsize = 1.5) {
  #' Paint every marker separately and save it in a folder
  #' 
  #' @param obj Seurat Object to paint with
  #' @param markers the vector markers to paint
  #' @param filepath the filepath to save the paintings
  
  obj@active.assay <- "RNA"
  for (i in 1:length(markers)) {
    gene = markers[i]
    png(filename = paste0(filepath, gene, ".png"), width = mywidth, height = myheight, unit="px")
    p = FeaturePlot(obj, features = c(gene), split.by = "sample", reduction = "umap", pt.size = myptsize, label=TRUE, order = TRUE)
    print(p)
    dev.off()
  } # end for
}

printVectorAsNewVector <- function(vect) {
  str = "newVect = c("
  for (i in 1:length(vect)) {
    element = vect[i]
    if (i == length(vect)) { # if it's the last element, don't put a comma
      str = paste0(str, '"', element, '")') 
    } else {
      str = paste0(str, '"', element, '", ') 
    }
  }
  str = paste0(str, "")
  cat(str)
  return(str)
}

bestCombo <- function(obj, targetCluster) {
  # Find the best combo of genes that are exclusive to that cluster
  Idents(obj) <- "seurat_clusters"
  targetCells <- WhichCells(obj, idents = targetCluster)
  targetGenes <- rownames(obj@assays$RNA@counts)[which(rowSums(obj@assays$RNA@counts[,targetCells]) > 100)]
  # maxTarget <- 0
  # minElse   <- ncol(obj)
  geneTargetCells <- list()
  geneOtherCells  <- list()
  for (i in targetGenes) {
    thisSum <- obj@assays$RNA@counts[i,]
    geneTargetCells[[i]] <- names(thisSum)[which(thisSum != 0 & names(thisSum) %in% targetCells)]
    geneOtherCells[[i]]  <- names(thisSum)[which(thisSum != 0 & ! names(thisSum) %in% targetCells)]
  }
  bestScore <- -ncol(obj)
  bestGenes <- c()
  n <- 0
  for (i in targetGenes) {
      if (n %% 1000 == 0) {
        print(n)
        print(bestScore)
      }
    for (j in targetGenes) {
      # print(paste(i,j))
      thisTarget <- length(geneTargetCells[[i]][which(geneTargetCells[[i]] %in% geneTargetCells[[j]])])
      thisOther  <- length(geneOtherCells[[i]][which(geneOtherCells[[i]] %in% geneOtherCells[[j]])])
      score <- thisTarget - thisOther
      if (score > bestScore) {
        bestScore <- score
        bestGenes <- c(i, j)
      }
    }
    n = n+1
  }
  # n <- 0
  # for (i in targetGenes) {
  #   if (n %% 1000 == 0) {
  #     print(n)
  #   }
  #   k <- 0
  #   for (j in targetGenes) {
  #     # print(k)
  #     thisSum <- colSums(obj@assays$RNA@counts[c(i,j),])
  #     thisTarget <- length(thisSum[which(thisSum != 0 & names(thisSum) %in% targetCells)])
  #     thisElse <- length(thisSum[which(thisSum != 0 & ! thisSum %in% targetCells)])
  #     if (thisTarget/thisElse > bestRatio) {
  #       # minElse <- thisElse
  #       # maxTarget <- thisTarget
  #       bestRatio <- thisTarget/thisElse
  #       bestGenes <- c(i, j)
  #     }
  #     k = k + 1
  #   }
  #   n = n + 1
  # }
}

myTotalTrans <- function(obj, slot="data", assay="RNA", features = NULL, cells = NULL) {
  result <- data.frame()
  for (ident in levels(Idents(obj))) {
    print(paste("Averaging Expression for", ident))
    if (identical(features, NULL)) {
      features = rownames(obj)
    }
    
    this_cells <- WhichCells(obj, idents = ident)
    if (! identical(cells, NULL)) {
      this_cells <- this_cells[which(this_cells %in% cells)]
    }
    
    if (length(features) == 1) {
      newRow <- setNames(as.data.frame(sum(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells])), features)
      names(newRow) <- features
      result <- rbind(result, newRow)
    } else {
      result <- rbind(result, setNames(t(as.data.frame(rowSums(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells]))), features)) 
    }
    rownames(result)[length(rownames(result))] <- ident
  }
  result <- as.data.frame(t(result))
  return(result)
}

myAverageExpression <- function(obj, slot="data", assay="RNA", features = NULL, cells=NULL) {
  result <- data.frame()
  for (ident in levels(Idents(obj))) {
    print(paste("Averaging Expression for", ident))
    if (identical(features, NULL)) {
      features = rownames(obj)
    }
    
    this_cells <- WhichCells(obj, idents = ident)
    if (! identical(cells, NULL)) {
      this_cells <- this_cells[which(this_cells %in% cells)]
    }
    
    if (length(features) == 1) {
      newRow <- setNames(as.data.frame(sum(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells])/length(this_cells)), features)
      names(newRow) <- features
      result <- rbind(result, newRow)
    } else if (length(this_cells) == 1) {
      newRow <- setNames(t(as.data.frame(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells])), features)
      names(newRow) <- features
      result = rbind(result, newRow)
    } else {
      newRow = setNames(t(as.data.frame(rowSums(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells]))/length(this_cells)), features)
      names(newRow) = features
      result <- rbind(result, newRow) 
    }
    rownames(result)[length(rownames(result))] <- ident
  }
  result <- as.data.frame(t(result))
  return(result)
}

convertHgncDataFrameToMzebra = function(df, gene_column, gene_names, na.rm = F, return_vect = F) {
  #' Convert a DataFram of HGNC gene names to ENS cichlid gene names
  #' 
  #'@param df dataframe of hgnc gene names
  #'@param gene_column column that has the hgnc gene names
  #'@param gene_names valid ENS cichlid gene names
  #'@param na.rm Remove NA ENS cichlid gene names
  #'@param return_vect return vector of ENS cichlid gene names instead of dataframe
  #'@return df data frame with a new column for ENS cichlid gene names called "mz" (or optionally return a vector)
  
  # Extract Input Genes
  hgnc_names = df[,gene_column]
  hgnc_names_unique = unique(hgnc_names)
  
  # Mart Objects
  mzebra = useEnsembl("ensembl", mirror = "uswest", dataset = "mzebra_gene_ensembl")
  human =  useEnsembl("ensembl", mirror = "uswest", dataset = "hsapiens_gene_ensembl")
  
  # Make a converter that includes the ENS names and Zebrafish names
  ens_converter = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = hgnc_names_unique , mart = human, attributesL = c("ensembl_gene_id", "zfin_id_symbol", "external_gene_name"), martL = mzebra, uniqueRows=T)
  ens_converter$Gene.name[which(ens_converter$Gene.name == "")] = NA
  
  # Make a converter that includes the ENS names, Zebrafish names, AND the input gene name to lower case.
  complete_converter = data.frame(hgnc = hgnc_names_unique)
  complete_converter$mz_ens = ens_converter$Gene.stable.ID[match(complete_converter$hgnc, ens_converter$HGNC.symbol)]
  complete_converter$mz_ens_name = ens_converter$Gene.name[match(complete_converter$hgnc, ens_converter$HGNC.symbol)]
  complete_converter$mz_zfin = ens_converter$ZFIN.symbol[match(complete_converter$hgnc, ens_converter$HGNC.symbol)]
  complete_converter$mz_lower = tolower(hgnc_names_unique)
  # complete_converter$mz_lower[which(! complete_converter$mz_lower %in% gene_names)] = NA
  complete_converter$mz_upper = hgnc_names_unique
  # complete_converter$mz_upper[which(! complete_converter$mz_upper %in% gene_names)] = NA
  complete_converter$mz_lower_many = paste(tolower(hgnc_names_unique), "(1 of many)")
  # complete_converter$mz_lower_many[which(! complete_converter$mz_lower_many %in% gene_names)] = NA
  complete_converter$mz_upper_many = paste(toupper(hgnc_names_unique), "(1 of many)")
  # complete_converter$mz_upper_many[which(! complete_converter$mz_upper_many %in% gene_names)] = NA
  
  # Pick which one of the possible cichlid gene names to use
  pickMZ = function(x) {
    good_mz = complete_converter[x, which(complete_converter[x,] %in% gene_names)[1]]
    if (length(good_mz) == 0)
      good_mz = NA
    return(good_mz)
  }
  complete_converter$mz = sapply(1:nrow(complete_converter), pickMZ)
  # complete_converter$mz = unlist(sapply(1:nrow(complete_converter), function(x) complete_converter[x, which(complete_converter[x,] %in% gene_names)[1]]))
  complete_converter$mz[which(complete_converter$mz == "NULL")] = NA
  
  # Convert Input Gene Names
  df$mz = complete_converter$mz[match(hgnc_names, complete_converter$hgnc)]
  if (na.rm)
    df = df[which(! is.na(df$mz) ),]
  if (return_vect)
    df = df$mz
  
  return(df)
}

convertToMouseObj <- function(obj) {
  # Converts a Mzebra Seurat object to a Mouse Seurat Object.
  
  # Convert a Seurat object to have all HGNC gene names
  print("Converting Genes Names...")
  genes <- rownames(obj@assays$RNA@counts)

  mouse  = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  
  # DF to convert from org to HGNC
  all_hgnc <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  
  # Initialize New Matricies
  print("Creating New Matrices...")
  new_counts_matrix <- as(obj@assays$RNA@counts, "sparseMatrix") 
  new_data_matrix   <- as(obj@assays$RNA@data, "sparseMatrix") 
  
  not_i <- 0
  multiple_hgnc <- 0
  bad_multiple_hgnc <- 0
  for (gene in genes) {
    if (gene %in% all_hgnc[,1]) {
      hgnc_gene <- all_hgnc[which(all_hgnc[,1] == gene),2]
      if (length(hgnc_gene) > 1) {
        multiple_hgnc <- multiple_hgnc + 1
        upper_hgnc_gene <- hgnc_gene[which(startsWith(tolower(hgnc_gene), tolower(gene)))]
        if (length(upper_hgnc_gene) == 1) {
          hgnc_gene <- upper_hgnc_gene
        } else {
          bad_multiple_hgnc <- bad_multiple_hgnc + 1
          hgnc_gene <- hgnc_gene[1]
        } # end bad multiple
      } # end multiple
      
      rownames(new_counts_matrix)[which(rownames(new_counts_matrix) == gene)] <- hgnc_gene
      rownames(new_data_matrix)[which(rownames(new_data_matrix) == gene)] <- hgnc_gene
    } else {
      not_i <- not_i + 1
    } # end gene not an hgnc gene
  } # end gene for
  print(paste("-Number of Genes not converted to HGNC:", not_i))
  print(paste("-Number of Genes with multple HGNC:", multiple_hgnc))
  print(paste("-Number of Genes with multple HGNC and non-ideal circumstance:", bad_multiple_hgnc))
  
  # Remove mzebra rows
  print("Removing old org rows...")
  ind <- which(rownames(new_counts_matrix) %in% all_hgnc[,2])
  new_counts_matrix <- new_counts_matrix[ind,]
  new_data_matrix   <- new_data_matrix[ind,]
  
  # Merge the duplicated rows
  print("Removing duplicated HGNC rows...")
  ptm <- proc.time()
  dup_genes <- unique(rownames(new_counts_matrix)[which(duplicated(rownames(new_counts_matrix)))])
  remove_dup_ind <- c()
  keep_dup_ind <- c()
  j <- 0
  dup_matrix_short_counts <- new_counts_matrix[which(rownames(new_counts_matrix) %in% dup_genes),]
  dup_matrix_short_data   <- new_data_matrix[which(rownames(new_data_matrix) %in% dup_genes),]
  keep_dup_matrix_counts  <- matrix(0L, nrow = length(dup_genes), ncol = ncol(obj@assays$RNA@counts))
  keep_dup_matrix_data    <- matrix(0L, nrow = length(dup_genes), ncol = ncol(obj@assays$RNA@data))
  rownames(keep_dup_matrix_counts) <- dup_genes
  rownames(keep_dup_matrix_data)   <- dup_genes
  for (gene in dup_genes) {
    ind_keep <- which(rownames(dup_matrix_short_counts) == gene)
    keep_dup_matrix_counts[j,] <- colSums(dup_matrix_short_counts[ind_keep,])
    keep_dup_matrix_data[j,]   <- colSums(dup_matrix_short_data[ind_keep,])
    j <- j + 1
  }
  # Delete all the duplicated rows at once
  print(paste("-Number of rows before merging:", nrow(new_counts_matrix)))
  print("-Actually doing the removal now I swear")
  remove_dup_ind <- which(rownames(new_counts_matrix) %in% dup_genes)
  new_counts_matrix <- new_counts_matrix[-remove_dup_ind,]
  new_data_matrix   <- new_data_matrix[-remove_dup_ind,]
  # nrow(new_counts_matrix)
  # Add all the new data at once
  print("-Combining the duplicated rows")
  new_counts_matrix <- rbind(new_counts_matrix, keep_dup_matrix_counts)
  new_data_matrix   <- rbind(new_data_matrix, keep_dup_matrix_data)
  print(paste("-Number of rows after merging:", nrow(new_counts_matrix)))
  proc.time() - ptm
  # new_counts_matrix[keep_dup_ind,] <- keep_dup_matrix_counts
  # new_data_matrix[keep_dup_ind,]   <- keep_dup_matrix_data
  # Delete temporary files
  rm(dup_matrix_short_counts)
  rm(dup_matrix_short_data)
  rm(keep_dup_matrix_counts)
  rm(keep_dup_matrix_data)
  
  print("Creating New Seurat Object...")
  obj_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj@project.name)
  obj_2 <- SetAssayData(object = obj_2, slot = 'data', new.data = new_data_matrix)
  # obj_2$seurat_clusters <- obj$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj@meta.data)) {
    obj_2@meta.data[col] <- obj@meta.data[col]
  }
  
  return(obj_2)
}

convertMzebraGeneListToMouse <- function(gene_list) {
  mzebra = useEnsembl("ensembl", mirror = "useast", dataset = "mzebra_gene_ensembl")
  mouse  = useEnsembl("ensembl", mirror = "useast", dataset = "mmusculus_gene_ensembl")
  
  # ensembl_genes <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = gene_list , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  # zfin_genes    <- getLDS(attributes = c("zfin_id_symbol"), filters = "zfin_id_symbol", values = gene_list , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  # hgnc_genes    <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = gene_list , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  # 
  # all_hgnc <- rbind(ensembl_genes, setNames(zfin_genes, names(ensembl_genes)), setNames(hgnc_genes, names(ensembl_genes)))
  # all_hgnc = all_hgnc[!duplicated(all_hgnc[,1]),]
  # all_mouse <- all_hgnc
  
  # DF to convert from org to HGNC
  all_mouse <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = gene_list, mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  
  mouse_genes <- unique(all_mouse[,2])
  print(paste0(length(mouse_genes)/length(gene_list) * 100, "% Genes Converted (", length(mouse_genes), "/", length(gene_list), ")"))
  return(all_mouse)
}

convertMzebraDFToMouse <- function(df, gene_column) {
  bad_genes <- df[,gene_column]
  bad_genes <- unique(bad_genes)
  gene_list <- bad_genes

  all_mouse <- convertMzebraGeneListToMouse(bad_genes)
  
  multiple_hgnc <- 0
  bad_multiple_hgnc <- 0
  converter <- data.frame()
  for (gene in bad_genes) {
    mouse_gene <- all_mouse[which(all_mouse[,1] == gene),2]
    new_mouse_gene <- mouse_gene[1]
    if (length(mouse_gene) > 1) {
      multiple_hgnc <- multiple_hgnc + 1
      upper_mouse_gene <- mouse_gene[which(startsWith(tolower(gene), tolower(mouse_gene)))]
      if (length(upper_mouse_gene) == 1) {
        new_mouse_gene <- upper_mouse_gene
      } else {
        bad_multiple_hgnc <- bad_multiple_hgnc + 1
        new_mouse_gene <- mouse_gene[1]
      } # end bad multiple
    } # end multiple
    converter <- rbind(converter, setNames(data.frame(t(c(gene, new_mouse_gene))), c("mzebra", "mouse")) )
  }
  print(paste0("Muliplte Mouse: ", multiple_hgnc))
  print(paste0("Bad Multiple Mouse: ", bad_multiple_hgnc))
  
  df[,gene_column] <- converter[match(df[,gene_column], converter[,1]),2]
  df <- df[which(! is.na(df[,gene_column])),]
  
  return(df)
}

keepCommonGenesObj <- function(obj_a, obj_b) {
  # Finds the common genes between two Seurat objects and
  # makes two new Seurat objects with just those genes.
  genes_a <- rownames(obj_a@assays$RNA@counts)
  genes_b <- rownames(obj_b@assays$RNA@counts)
  common <- genes_a[which(genes_a %in% genes_b)]
  print(paste("Found", length(common), "genes in common -", length(common)/length(genes_a) * 100, "of obj_a genes", length(common)/length(genes_b) * 100, "of obj_b genes" ))
  
  ## First Seurat Object
  # Initialize New Matricies
  print("Creating New Matrices (for First Seurat Object)...")
  new_counts_matrix <- as(obj_a@assays$RNA@counts, "sparseMatrix") 
  new_data_matrix   <- as(obj_a@assays$RNA@data, "sparseMatrix")
  
  # Removing Non-Overlapping Genes
  print("Removing Non-Overlapping Genes (for the First Seurat Object)...")
  all_ind_keep <- c()
  all_ind <- 1:nrow(new_counts_matrix)
  for (gene in common) {
    ind_keep <- which(rownames(new_counts_matrix) == gene)
    all_ind_keep <- c(all_ind_keep, ind_keep)
  }
  all_ind_remove <- all_ind[which(! all_ind %in% all_ind_keep)]
  new_counts_matrix <- new_counts_matrix[-all_ind_remove,]
  new_data_matrix   <- new_data_matrix[-all_ind_remove,]
  
  print("Creating New Seurat Object (for the First Seurat Object)...")
  obj_a_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj_a@project.name)
  obj_a_2 <- SetAssayData(object = obj_a_2, slot = 'data', new.data = new_data_matrix)
  # obj_a_2$seurat_clusters <- obj_a$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj_a@meta.data)) {
    obj_a_2@meta.data[col] <- obj_a@meta.data[col]
  }
  
  ## Second Seurat Object
  # Initialize New Matricies
  print("Creating New Matrices (for Second Seurat Object)...")
  new_counts_matrix <- as(obj_b@assays$RNA@counts, "sparseMatrix") 
  new_data_matrix   <- as(obj_b@assays$RNA@data, "sparseMatrix")
  
  # Removing Non-Overlapping Genes
  print("Removing Non-Overlapping Genes (for the Second Seurat Object)...")
  all_ind_keep <- c()
  all_ind <- 1:nrow(new_counts_matrix)
  for (gene in common) {
    ind_keep <- which(rownames(new_counts_matrix) == gene)
    all_ind_keep <- c(all_ind_keep, ind_keep)
  }
  all_ind_remove <- all_ind[which(! all_ind %in% all_ind_keep)]
  new_counts_matrix <- new_counts_matrix[-all_ind_remove,]
  new_data_matrix   <- new_data_matrix[-all_ind_remove,]
  
  print("Creating New Seurat Object (for the Second Seurat Object)...")
  obj_b_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj_b@project.name)
  obj_b_2 <- SetAssayData(object = obj_b_2, slot = 'data', new.data = new_data_matrix)
  # obj_b_2$seurat_clusters <- obj_b$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj_b@meta.data)) {
    obj_b_2@meta.data[col] <- obj_b@meta.data[col]
  }
  
  return(list(obj_a_2, obj_b_2))
}

convertToHgncObj <- function(obj, organism) {
  # Convert a Seurat object to have all HGNC gene names
  print("Converting Genes Names...")
  genes <- rownames(obj@assays$RNA@counts)
  if (organism == "mouse") {
    dataset_name <- "mmusculus_gene_ensembl"
  } else if (organism == "mzebra") {
    dataset_name <- "mzebra_gene_ensembl"
  } else {
    stop("Organism not recognized. Options are mouse or mzebra.")
  }
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  org = useMart(biomart="ensembl", dataset=dataset_name)
  
  # DF to convert from org to HGNC
  all_hgnc <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = org, attributesL = c("external_gene_name"), martL = human, uniqueRows=T)
  
  # Initialize New Matricies
  print("Creating New Matrices...")
  new_counts_matrix <- as(obj@assays$RNA@counts, "sparseMatrix") 
  new_data_matrix   <- as(obj@assays$RNA@data, "sparseMatrix") 
  
  not_i <- 0
  multiple_hgnc <- 0
  bad_multiple_hgnc <- 0
  for (gene in genes) {
    if (gene %in% all_hgnc[,1]) {
      hgnc_gene <- all_hgnc[which(all_hgnc[,1] == gene),2]
      if (length(hgnc_gene) > 1) {
        multiple_hgnc <- multiple_hgnc + 1
        upper_hgnc_gene <- hgnc_gene[which(startsWith(tolower(hgnc_gene), tolower(gene)))]
        if (length(upper_hgnc_gene) == 1) {
          hgnc_gene <- upper_hgnc_gene
        } else {
          bad_multiple_hgnc <- bad_multiple_hgnc + 1
          hgnc_gene <- hgnc_gene[1]
        } # end bad multiple
      } # end multiple
      
      rownames(new_counts_matrix)[which(rownames(new_counts_matrix) == gene)] <- hgnc_gene
      rownames(new_data_matrix)[which(rownames(new_data_matrix) == gene)] <- hgnc_gene
    } else {
      not_i <- not_i + 1
    } # end gene not an hgnc gene
  } # end gene for
  print(paste("-Number of Genes not converted to HGNC:", not_i))
  print(paste("-Number of Genes with multple HGNC:", multiple_hgnc))
  print(paste("-Number of Genes with multple HGNC and non-ideal circumstance:", bad_multiple_hgnc))
  
  # Remove mzebra rows
  print("Removing old org rows...")
  ind <- which(rownames(new_counts_matrix) %in% all_hgnc[,2])
  new_counts_matrix <- new_counts_matrix[ind,]
  new_data_matrix   <- new_data_matrix[ind,]
  
  # Merge the duplicated rows
  print("Removing duplicated HGNC rows...")
  ptm <- proc.time()
  dup_genes <- unique(rownames(new_counts_matrix)[which(duplicated(rownames(new_counts_matrix)))])
  remove_dup_ind <- c()
  keep_dup_ind <- c()
  j <- 0
  dup_matrix_short_counts <- new_counts_matrix[which(rownames(new_counts_matrix) %in% dup_genes),]
  dup_matrix_short_data   <- new_data_matrix[which(rownames(new_data_matrix) %in% dup_genes),]
  keep_dup_matrix_counts  <- matrix(0L, nrow = length(dup_genes), ncol = ncol(obj@assays$RNA@counts))
  keep_dup_matrix_data    <- matrix(0L, nrow = length(dup_genes), ncol = ncol(obj@assays$RNA@data))
  rownames(keep_dup_matrix_counts) <- dup_genes
  rownames(keep_dup_matrix_data)   <- dup_genes
  for (gene in dup_genes) {
    ind_keep <- which(rownames(dup_matrix_short_counts) == gene)
    keep_dup_matrix_counts[j,] <- colSums(dup_matrix_short_counts[ind_keep,])
    keep_dup_matrix_data[j,]   <- colSums(dup_matrix_short_data[ind_keep,])
    j <- j + 1
  }
  # Delete all the duplicated rows at once
  print(paste("-Number of rows before merging:", nrow(new_counts_matrix)))
  print("-Actually doing the removal now I swear")
  remove_dup_ind <- which(rownames(new_counts_matrix) %in% dup_genes)
  new_counts_matrix <- new_counts_matrix[-remove_dup_ind,]
  new_data_matrix   <- new_data_matrix[-remove_dup_ind,]
  # nrow(new_counts_matrix)
  # Add all the new data at once
  print("-Combining the duplicated rows")
  new_counts_matrix <- rbind(new_counts_matrix, keep_dup_matrix_counts)
  new_data_matrix   <- rbind(new_data_matrix, keep_dup_matrix_data)
  print(paste("-Number of rows after merging:", nrow(new_counts_matrix)))
  tot_before = sum(rowSums(obj@assays$RNA@counts))
  tot_after = sum(rowSums(new_counts_matrix))
  print(paste0("-Number of counts lost: ", tot_before-tot_after, " (", format(round( (tot_before-tot_after)/tot_before * 100 , 2), nsmall = 2), "%)"))
  proc.time() - ptm
  # new_counts_matrix[keep_dup_ind,] <- keep_dup_matrix_counts
  # new_data_matrix[keep_dup_ind,]   <- keep_dup_matrix_data
  # Delete temporary files
  rm(dup_matrix_short_counts)
  rm(dup_matrix_short_data)
  rm(keep_dup_matrix_counts)
  rm(keep_dup_matrix_data)
  
  print("Creating New Seurat Object...")
  obj_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj@project.name)
  obj_2 <- SetAssayData(object = obj_2, slot = 'data', new.data = new_data_matrix)
  # obj_2$seurat_clusters <- obj$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj@meta.data)) {
    obj_2@meta.data[col] <- obj@meta.data[col]
  }
  
  return(obj_2)
}

geneCap <- function(gene, gene_names) {
  # Gene the gene name in the right format
  gene_lower <- tolower(gene)
  gene_upper <- toupper(gene)
  gene_title <- str_to_title(gene)
  error <- FALSE
  if (gene_lower %in% gene_names) {
    gene <- gene_lower
  } else if (gene_upper %in% gene_names) {
    gene <- gene_upper
  } else if (gene_title %in% gene_names) {
    gene <- gene_title
  } else {
    error <- TRUE
  }
  
  return(c(gene, error))
}

validGenes <- function(genes, gene_names) {
  valid_genes <- c()
  for (gene in genes) {
    result <- geneCap(gene, gene_names)
    gene <- result[1]
    error <- as.logical(result[2])
    if (! error) {
      valid_genes <- c(valid_genes, gene)
    }
  } # end gene for
  valid_genes <- unique(valid_genes)
  return(valid_genes)
} # end validGenes function


convertMouseDataFrameToHgnc = function(mouse_df, gene_column) {
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse  = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  m_genes = mouse_df[,gene_column]
  
  hgnc_df = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = m_genes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  colnames(hgnc_df) = c("mouse", "hgnc")
  hgnc_df = hgnc_df[which(! is.na(hgnc_df$hgnc)),]
  
  converter = hgnc_df
  mouse_df[,gene_column] <- converter[match(mouse_df[,gene_column], converter[,1]),2]
  mouse_df <- mouse_df[which(! is.na(mouse_df[,gene_column])),]
  
  return(mouse_df)
}

convertToHgnc <- function(genes) {
  # human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  mzebra = useEnsembl("ensembl", mirror = "uswest", dataset = "mzebra_gene_ensembl")
  human =  useEnsembl("ensembl", mirror = "uswest", dataset = "hsapiens_gene_ensembl")
  
  ensembl_genes <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  zfin_genes    <- getLDS(attributes = c("zfin_id_symbol"), filters = "zfin_id_symbol", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  hgnc_genes    <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  all_hgnc <- rbind(ensembl_genes, setNames(zfin_genes, names(ensembl_genes)), setNames(hgnc_genes, names(ensembl_genes)))
  all_hgnc = all_hgnc[!duplicated(all_hgnc[,1]),]
  colnames(all_hgnc) <- c("mzebra", "hgnc")
  
  # all_hgnc <- unique(c(ensembl_genes[,2], zfin_genes[,2], hgnc_genes[,2]))
  return(all_hgnc)
}

hgncMzebra <- function(genes, gene_names, onPACE=F) {
  # genes: genes to find hgnc conversions for
  # gene_names: list of ALL mzebra gene names
  # onPACE: is the script being run on PACE (default = FALSE)
  if (onPACE) {
    pat <- read.table("~/scratch/m_zebra_ref/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE)
  } else {
    # pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE)
    pat = read.table("~/research/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE)
  }
  valid_genes <- validGenes(genes, gene_names)
  all_hgnc <- convertToHgnc(gene_names)
  ind_bool <- match(genes,pat$V2)
  ind <- ind_bool[! is.na(ind_bool)]
  found_names <- as.vector(pat$V7[ind])
  found_names <- found_names[!is.na(found_names)]
  found_names_hgnc <- as.vector(pat$V8[ind])
  # found_names_hgnc <- found_names_hgnc[!is.na(found_names_hgnc)]
  
  df1 <- setNames(as.data.frame(found_names_hgnc), c("hgnc"))
  found_names_hgnc = setNames(data.frame(genes[which(! is.na(ind_bool))], found_names_hgnc), c("mzebra", "hgnc"))
  found_names_hgnc = found_names_hgnc[which(! is.na(found_names_hgnc$hgnc)),]
  # found_names_hgnc <- inner_join(df1, all_hgnc, by = "hgnc")
  
  pseudo_hgnc <- toupper(genes)
  df1 <- setNames(as.data.frame(pseudo_hgnc), c("hgnc"))
  found_mzebra <- inner_join(df1, all_hgnc, by = "hgnc")
  
  found_mzebra <- found_mzebra[,2:1]
  # found_names_hgnc <- found_names_hgnc[,2:1]
  good_df <- rbind(all_hgnc, setNames(found_names, names(all_hgnc)), setNames(found_mzebra, names(all_hgnc)), setNames(found_names_hgnc, names(all_hgnc)))
  good_df <- unique(good_df)
  return(good_df)
}

hgncMouse <- function(genes) {
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse  = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  
  ensembl_genes <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(ensembl_genes)
}

hgncMzebraInPlace <- function(df, gene_column, gene_names, onPACE=F, rm_na = T) {
  bad_genes <- df[,gene_column]
  bad_genes <- unique(bad_genes)
  converter <- hgncMzebra(bad_genes, gene_names, onPACE)
  df[,gene_column] <- converter[match(df[,gene_column], converter[,1]),2]
  if (rm_na) {
    df <- df[which(! is.na(df[,gene_column])),]
  }
  
  return(df)
}

hgncMouseInPlace <- function(df, gene_column) {
  bad_genes <- df[,gene_column]
  bad_genes <- unique(bad_genes)
  converter <- hgncMouse(bad_genes)
  df[,gene_column] <- converter[match(df[,gene_column], converter[,1]),2]
  df <- df[which(! is.na(df[,gene_column])),]
  
  return(df)
}

# THIS IS THE UPDATED/GOOD Function 02/21/2020
hgncGood <- function(genes, gene_names, as_df = F) {
  pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE)
  valid_genes <- validGenes(genes, gene_names)
  all_hgnc <- convertToHgnc(gene_names)
  ind <- match(genes,pat$V2)
  ind <- ind[! is.na(ind)]
  found_names <- as.vector(pat$V7[ind])
  found_names <- found_names[!is.na(found_names)]
  found_names_hgnc <- as.vector(pat$V8[ind])
  found_names_hgnc <- found_names_hgnc[!is.na(found_names_hgnc)]
  
  df1 <- setNames(as.data.frame(found_names_hgnc), c("hgnc"))
  found_names_hgnc <- inner_join(df1, all_hgnc, by = "hgnc")
  if (! as_df) {
    found_names_hgnc <- found_names_hgnc$mzebra
  }
  # ind_found_hgnc <- match(found_names_hgnc, all_hgnc$hgnc)
  # ind_found_hgnc <- ind_found_hgnc[! is.na(ind_found_hgnc)]
  # found_names_hgnc <- as.vector(all_hgnc[ind_found_hgnc,1])
  
  pseudo_hgnc <- toupper(genes)
  df1 <- setNames(as.data.frame(pseudo_hgnc), c("hgnc"))
  found_mzebra <- inner_join(df1, all_hgnc, by = "hgnc")
  if (! as_df) {
    found_mzebra = found_mzebra$mzebra
  }
  # ind_hgnc <- match(pseudo_hgnc, all_hgnc$hgnc)
  # ind_hgnc <- ind_hgnc[! is.na(ind_hgnc)]
  # found_mzebra <- as.vector(all_hgnc[ind_hgnc,1])
  
  if (! as_df) {
    good <- c(valid_genes, found_names, found_names_hgnc, found_mzebra)
    good <- good[which(good != "")]
    good <- validGenes(good, gene_names)
    good <- unique(good)
    good <- sort(good)
  } else {
    valid_genes_df = data.frame(hgnc = toupper(valid_genes), mzebra = valid_genes)
    good = rbind(found_names_hgnc, found_mzebra, valid_genes_df)
    good = good[which(good$mzebra != ""),]
    good = good[which(good$mzebra %in% validGenes(good$mzebra, gene_names)),]
    good = good[which(! duplicated(good$mzebra)),]
  }
  
  return(good)
}

# This function is outdated
onlyGood <- function(genes, gene_names) {
  pat <- read.table("C:/Users/miles/Downloads/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE)
  valid_genes <- validGenes(genes, gene_names)
  ind <- match(genes,pat$V2)
  ind <- ind[! is.na(ind)]
  found_names <- as.vector(pat$V7[ind])
  found_names <- found_names[!is.na(found_names)]
  found_names_hgnc <- as.vector(pat$V8[ind])
  found_names_hgnc <- found_names_hgnc[!is.na(found_names_hgnc)]
  good <- c(valid_genes, found_names, found_names_hgnc)
  good <- unique(good)
  good <- sort(good)
  good <- good[which(good != "")]
  return(good)
}

goodInPlace <- function(data, gene_column, gene_names) {
  # Only keep the rows in the dataframe that have gene names that are in our data
  pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE)
  keep_rows <- c()
  for (i in 1:nrow(data)) {
    gene <- as.character(data[i,gene_column])
    result <- geneCap(gene, gene_names)
    gene <- result[1]
    error <- as.logical(result[2])
    ens2_gene <- as.character(pat$V6[which(gene == pat$V2)])[1]
    if (!error) {
      data[i,gene_column] <- gene
      keep_rows <- c(keep_rows, i)
    }
    if (!identical(ens2_gene, character(0)) && ens2_gene %in% gene_names) {
      data[i,gene_column] <- ens2_gene
      keep_rows <- c(keep_rows, i)
    }
  }
  new_data <- data[keep_rows,]
  return(new_data)
}

ncAddAllGeneInfo = function(df, gene_column, onPACE = F) {
  #' Adds all the gene info for a gene to the dataframe, including:
  #' Ensembl name, human ortholog, biotype, description, human ortholog description
  #' 
  #' @param df input dataframe with a column with MZ NCBI genes
  #' @param gene_column column number that contains MZ NCBI genes
  #' @param onPACE is the code being run on PACE?
  
  if (onPACE) {
    gene_info = read.table("~/scratch/m_zebra_ref/gene_info.txt", sep="\t", header = T, stringsAsFactors = F)
  } else {
    gene_info = read.table("C:/Users/miles/Downloads/all_research/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
  }
  
  genes = df[,gene_column]
  df$ens = gene_info$ens[match(genes, gene_info$mzebra)]
  df$human = gene_info$human[match(genes, gene_info$mzebra)]
  df$biotype = gene_info$biotype[match(genes, gene_info$mzebra)]
  df$mz_description = gene_info$mzebra_description[match(genes, gene_info$mzebra)]
  df$human_description = gene_info$human_description[match(genes, gene_info$mzebra)]
  
  return(df)
}

ncEnsGeneConverter = function(genes, isNc = T, onPACE = F, rm_na = T) {
  #' Converts from NCBI mzebra genes to ENS mzebra genes and vice versa.
  #' 
  #' @param genes input vector of genes
  #' @param isNc is the input list of genes in NCBI format?
  #' @param onPACE is the code being run on PACE?
  #' @param rm_na remove genes that weren't converted?
  #' @return convereted a vector of converted genes
  
  if (onPACE) {
    pat <- read.table("~/scratch/m_zebra_ref/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = F)
  } else {
    pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = F)
  }
  
  if (isNc) {
    converted = pat[match(genes, pat[,2]),7]
    converted[which( is.na(converted) & genes %in% pat[,7] )] = genes[which( is.na(converted) & genes %in% pat[,7] )]
  } else {
    converted = pat[match(genes, pat[,7]),2]
    converted[which( is.na(converted) & genes %in% pat[,2] )] = genes[which( is.na(converted) & genes %in% pat[,2] )]
  }
  
  if (rm_na)
    converted = converted[which( ! is.na(converted) )]
  
  return(converted)
}

ncAddDescription = function(df, gene_column, rm_na = F) {
  #' Given a dataframe, add a description to the NCBI mzebra genes
  #' 
  #' @param df input dataframe
  #' @param gene_column the column of the dataframe containing the mzebra genes
  #' @param rm_na remove rows where no description was found?
  #' @return df output dataframe with a description as the last column
  
  if (onPACE) {
    gene_info = read.table("~/scratch/m_zebra_ref/gene_info.txt", sep="\t", header = T, stringsAsFactors = F)
  } else {
    gene_info = read.table("C:/Users/miles/Downloads/all_research/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
  }
  
  genes = as.vector(df[,gene_column])
  
  df$description = gene_info$mzebra_description[match(genes, gene_info$mzebra)]
  if(rm_na)
    df = df[which(! is.na(df$description)),]
  
  return(df)
}

ncHgncMzebraInPlace = function(df, gene_column, rm_na = T, onPACE = F) {
  #' Given a datafrme, convert the NCBI mzebra genes in the specified column to NCBI.
  #' 
  #' @param df input dataframe
  #' @param gene_column the column of the dataframe containing the mzebra genes
  #' @param rm_na remove rows where no HGNC name was found for the mzebra gene?
  #' @param onPACE is the code being run on PACE?
  #' @return df a datafram with a hgnc column
  
  genes = as.vector(df[,gene_column])
  converter = ncMzToHgncVect(genes)
  df[,gene_column] = converter$human[match(df[,gene_column], converter$mzebra)]
  if (rm_na)
    df = df[which(! is.na(df[,gene_column])),]
  
  return(df)
}

# ncMzToHgncVectShortcut(mz, returnDF = T, onPACE=F, rm_na=F) {
  #' Convert NCBI mzebra gens to HGNC format. 
  #' But use file that I've already curated to save time.
  #' 
# }

ncMzToHgncVect = function(mz, returnDF = T, onPACE=F, rm_na = F) {
  #' Convert NCBI mzebra genes to HGNC format.
  #'
  #' @param mz input vector of genes in NCBI mzebra format
  #' @param returnDF return the output in the form of a dataframe? Else just a vector HGNC genes are returned
  #' @param onPACE is the code being run on PACE? The filepath for Patrick's file is different
  #' @param rm_na remove rows where no HGNC gene was found?
  #' @return all_hgnc a dataframe by default with 4 columns. First column is the mzebra gene,
  #' second column is the HGNC gene found in Patrick's file, the third column is the HGNC gene
  #' found by biomaRt, and the last column is my pick for the final HGNC name. Sometimes Patrick
  #' and biomaRt disagree, in those cases I keep Patrick's, unless the biomaRt gene is the 
  #' mzebra gene capitalized.
  
  if (onPACE) {
    pat <- read.table("~/scratch/m_zebra_ref/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = F)
  } else {
    pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = F)
  }
  
  ind_bool = match(mz,pat$V2)
  ind = ind_bool[! is.na(ind_bool)]
  pat_hgnc = data.frame(mzebra=as.vector(pat$V2[ind]), human_pat=as.vector(pat$V8[ind]))
  
  mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mart_hgnc = getLDS(attributes = c("entrezgene_accession"), filters = "entrezgene_accession", values = mz, mart = mzebra, attributesL = c("external_gene_name"), martL = human, uniqueRows=T)
  colnames(mart_hgnc) = c("mzebra", "human_mart")
  
  all_hgnc = full_join(pat_hgnc, mart_hgnc, by = "mzebra")
  all_hgnc$mzebra = as.vector(all_hgnc$mzebra)
  all_hgnc$human_pat = as.vector(all_hgnc$human_pat)
  all_hgnc$human_mart = as.vector(all_hgnc$human_mart)
  if (rm_na)
    all_hgnc[which(! is.na(all_hgnc$human_pat) | ! is.na(all_hgnc$human_mart)),]
  all_hgnc$human = all_hgnc$human_pat
  all_hgnc$human[which(! is.na(all_hgnc$human_mart)
                       & is.na(all_hgnc$human_pat))] = all_hgnc$human_mart[which(! is.na(all_hgnc$human_mart) & is.na(all_hgnc$human_pat))]
  disagree_mart_better = which(! is.na(all_hgnc$human_mart) 
                               & all_hgnc$human_mart != all_hgnc$human_pat
                               & (all_hgnc$human_mart == toupper(all_hgnc$mzebra) | all_hgnc$human_mart == toupper(substr(all_hgnc$mzebra, 1, nchar(all_hgnc$mzebra)-1))))
  all_hgnc$human[disagree_mart_better] = all_hgnc$human_mart[disagree_mart_better]
  if (! returnDF)
    all_hgnc = all_hgnc$human
  
  return(all_hgnc)
}

umd1To2a = function(genes, onPACE = F, as_df = F, rm_na = T) {
  #' Find the UMD2a list of genes from the UMD1 list
  #' 
  #' @param genes input UMD1 genes that will be converted
  #' @param onPACE is the code being run on PACE?
  #' @param as_df return the output as dataframe? With one column for UMD1 genes and one for the gene converted to UMD2a.
  #' @param rm_na remove UMD1 genes not converted
  #' @return umd2a vector of UMD2a genes
  
  if (onPACE) {
    pat <- read.table("~/scratch/m_zebra_ref/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = F)
  } else {
    pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = F)
  }
  
  genes = as.vector(genes)
  df = data.frame(umd1=genes, umd2a = NA)
  df$umd1 = as.vector(df$umd1)
  df$umd2a = as.vector(df$umd2a)
  
  df$umd2a[which(df$umd1 %in% pat$V2)] = df$umd1[which(df$umd1 %in% pat$V2)]
  df$umd2a[which(is.na(df$umd2a))] = pat$V2[match(df$umd1[which(is.na(df$umd2a))], pat$V1)]
  df$umd2a[which(is.na(df$umd2a))] = pat$V2[match(df$umd1[which(is.na(df$umd2a))], pat$V30)]
  df$umd2a = unlist(df$umd2a)
  
  pct = format(round(length(which(!is.na(df$umd2a)))/length(genes) * 100, 2), nsmall = 2)
  print(paste0( length(which(!is.na(df$umd2a))), " out of ", length(genes), " UMD1 genes were converted to UMD2a genes (", pct, "%)" ))
  if (rm_na) { df = df[which(!is.na(df$umd2a)),] }
  if (as_df) { return(df) }
  
  return(df$umd2a)
}

expressionDend = function(objs, my_slot="counts", imp_genes = NA) {
  # Purpose: Create dendrograms of expression of genes in clusters of mulitple objects
  # Input:
  #       obj: list of Seurat objects
  # Output: dendrograms
  
  # 1. Find Genes to Use in Dendrogam (using all genes would be too many).
  non_zero_genes = c()
  all_clusters = c()
  for (obj in objs) {
    gene_names = rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) != 0)]
    non_zero_genes = c(non_zero_genes, gene_names)
    all_clusters = c(all_clusters, paste(obj$project[[1]], sort(unique(as.vector((obj$seurat_clusters))))))
  }
  non_zero_table = table(non_zero_genes)
  non_zero_genes = names(non_zero_table)[which(non_zero_table == length(objs))]
  
  # imp_genes = c()
  # min_pct = 0.1
  # for (obj in objs) {
  #   print(obj$project[[1]])
  #   Idents(obj) = obj$seurat_clusters
  #   mat = obj@assays$RNA@counts
  #   mat[which(mat > 1)] = 1
  #   for (cluster in levels(obj$seurat_clusters)) {
  #     cells_cluster = WhichCells(object = obj, idents = cluster)
  #     n_cells_min = min_pct * length(cells_cluster)
  #     genes_pass = rownames(obj)[which(rowSums(as.matrix(mat[non_zero_genes,cells_cluster])) >= n_cells_min)]
  #     imp_genes = c(imp_genes, genes_pass)
  #   } # end cluster for
  # } # end obj for
  
  if (is.na(imp_genes)) {
    print("Genes Not Supplied. Using Non-Zero Genes to Use in Dendrogram.")
    imp_genes = non_zero_genes
  }
  imp_genes = unique(imp_genes)
  print(paste("Using", length(imp_genes), "in dendrogram"))
  
  # 2. Get Expression Data at Genes
  dend_mat = matrix(0L, nrow=length(imp_genes), ncol=length(all_clusters), dimnames = list(imp_genes, all_clusters))
  print(paste("Creating Dendrogram Matrix of size", nrow(dend_mat), "x", ncol(dend_mat)))
  for (obj in objs) {
    print(obj$project[[1]])
    # if ("annot" %in% colnames(obj@meta.data)) {
    #   clusters = unique(obj$annot)
    #   Idents(obj) = obj$annot
    # } else {
    #   clusters = levels(obj$seurat_clusters)
    #   Idents(obj) = obj$seurat_clusters
    # }
    Idents(obj) = obj$seurat_clusters
    clusters = sort(unique(as.vector(obj$seurat_clusters)))
    for (cluster in clusters) {
      cells_cluster = WhichCells(object = obj, idents = cluster)
      cluster_name = paste(obj$project[[1]], cluster)
      dend_mat[imp_genes,cluster_name] = rowMeans(obj@assays$RNA@data[imp_genes, cells_cluster])
    }
  }
  
  # dend_mat[which(dend_mat > 4)] = 4
  # dend_mat = log2(dend_mat)
  # dend_mat <- dend_mat[is.finite(rowSums(dend_mat)),]
  
  # Create the Plot
  png5_name = "test_dend.png"
  Colors=rev(brewer.pal(11,"Spectral"))
  # my_breaks = c(seq(-10, -2, length=100), seq(-2,-0.1,length=100), seq(-0.1, 0.1, length=100), seq(0.1,2,length=100), seq(2,10,length=100))
  my_palette1=colorRampPalette(Colors)(n = 299)
  my_palette <- colorRampPalette(c("#ff4b5c", "#FFFFFF", "#056674"))(n = 299)
  print("Plotting the dendrogram")
  par(mar=c(10, 4.1, 4.1, 2.1))
  png(png5_name, width = 50*length(all_clusters)+50, height = 50*length(all_clusters), unit = "px", res = 120)
  heatmap.2(dend_mat, scale = "none", dendrogram = "both", col = my_palette1, trace = "none", margins=c(10,5), srtCol=45, symm=F,symkey=F,symbreaks=F, main="Mean Expression")
  dev.off()
  print("Finished plotting")
  
  # genes_0 = c()
  # for (row in rownames(dend_mat)) {
  #   if (all(dend_mat[row,] == 0))
  #     genes_0 = c(genes_0, row)
  # }
  # length(genes_0)
  
  return(TRUE)
}

heatmapComparisonMulti = function(dfs, samples, filename, filepath, correction_factors, labels=F, xlab=T, tri=T, dendrogram=T) {
  # Input: list of dataframes that are output of Seurat FindAllMarkers
  #        Vector of samples or whatever you want to name those two dataframes
  #        Base file name for png output
  #        Filepath for png output
  #        org: a vector of organisms the samples belong to ("mzebra", "mouse", "human")
  clusters = list()
  num_clusters = list()
  all_logFC = c()
  all_genes = c()
  all_clusters = c()
  with_correction = ""
  for (i in 1:length(dfs)) {
    clusters[[i]] = unique(as.vector(dfs[[i]]$cluster))
    if (! any(is.na(as.numeric(clusters[[i]])))) {
      print("Can sort numerically")
      clusters[[i]] = sort(as.numeric(clusters[[i]]))
    }
    num_clusters[[i]] = length(clusters[[i]])
    all_logFC = c(all_logFC, dfs[[i]]$avg_logFC)
    all_clusters = c(all_clusters, paste(samples[[i]], clusters[[i]]))
  }
  
  # # Correct for gene conversion
  # correction_factor = list()
  # correction_factor[["mzebra"]] = 28622/16735  # all mzebra genes, all converted mzebra genes
  # correction_factor[["mouse"]]  = 23420/17140 
  # correction_factor[["human"]]  = 1
  
  # Now do Pairwise Comparison of each df's DEGs
  df = data.frame() # big df of all pairwise comparisons
  for (i in 1:length(dfs)) {
    for (i_clust in 1:num_clusters[[i]]) {
      # print(paste0("SAMPLE 1:", clusters[[i]][[i_clust]]))
      i_clust_df = dfs[[i]][which(dfs[[i]]$cluster == clusters[[i]][i_clust]),]
      i_clust_df = i_clust_df[!duplicated(i_clust_df$gene),]
      
      for (j in 1:length(dfs)) {
        for (j_clust in 1:num_clusters[[j]]) {
          # print(paste0("SAMPLE 2:", clusters[[j]][[j_clust]]))
          j_clust_df = dfs[[j]][which(dfs[[j]]$cluster == clusters[[j]][j_clust]),]
          j_clust_df = j_clust_df[!duplicated(j_clust_df$gene),]
          
          ovlp_genes = unique(j_clust_df$gene[which(j_clust_df$gene %in% i_clust_df$gene)])
          ovlp_genes = ovlp_genes[which(! is.na(ovlp_genes))]
          ovlp = length(ovlp_genes)
          j_clust_sign = sign(j_clust_df$avg_logFC[which(j_clust_df$gene %in% ovlp_genes)])
          i_clust_sign = sign(i_clust_df$avg_logFC[which(i_clust_df$gene %in% ovlp_genes)])
          ovlp_same_dir_genes = unique(ovlp_genes[which(j_clust_sign == i_clust_sign)])
          ovlp_same_dir = length(ovlp_same_dir_genes)
          all_genes = c(all_genes, ovlp_same_dir_genes)
          
          # Rename the clusters with their sample names to avoid confusion
          sample1_clust = paste0(samples[[i]], " ", clusters[[i]][i_clust])
          sample2_clust = paste0(samples[[j]], " ", clusters[[j]][j_clust])
          
          total_ovlp = 2*ovlp
          total_ovlp_same_dir = 2*ovlp_same_dir
          if ("correction_factor" %in% colnames(j_clust_df) && "correction_factor" %in% colnames(i_clust_df)) {
            total_ovlp = ovlp*i_clust_df$correction_factor[1] + ovlp*j_clust_df$correction_factor[1]
            total_ovlp_same_dir = ovlp_same_dir*i_clust_df$correction_factor[1] + ovlp_same_dir*j_clust_df$correction_factor[1]
            with_correction = "w/ Correction for Gene Conversion"
          }
          
          pct = (total_ovlp / (nrow(i_clust_df) + nrow(j_clust_df))) * 100
          pct_same_dir = (total_ovlp_same_dir / (nrow(i_clust_df) + nrow(j_clust_df))) * 100
          
          # print(sample1_clust)
          # print(sample2_clust)
          # print(pct_same_dir)
          # Check if pct is greater than 100
          if (sample1_clust != sample2_clust && pct_same_dir > 100) {
            print("Error pct ovlp > 100")
            print(sample1_clust)
            print(sample2_clust)
            print(ovlp_same_dir)
            print(nrow(i_clust_df))
            print(nrow(j_clust_df))
          }
          
          df <- rbind(df, t(c(sample1_clust, sample2_clust, ovlp, pct, ovlp_same_dir, pct_same_dir)))
        }
      }
      
    }
    print(paste("Finished Pairwise Comparisons for", samples[[i]]))
  } # end parwise comparison
  
  colnames(df) <- c("df1_cluster", "df2_cluster", "ovlp", "pct", "ovlp_same_dir", "pct_same_dir")
  df$ovlp = as.numeric(as.vector(df$ovlp))
  df$pct = as.numeric(as.vector(df$pct))
  df$ovlp_same_dir = as.numeric(as.vector(df$ovlp_same_dir))
  df$pct_same_dir = as.numeric(as.vector(df$pct_same_dir))
  
  # Extract lower triangle
  if (tri) {
    print("Only Keeping Lower Triangle")
    new_df = data.frame()
    clusters = unique(as.vector(df$df1_cluster))
    clusters = clusters[which(! is.na(clusters))]
    for (i in 1:length(clusters)) {
      i_clust = clusters[i]
      # if (i == 1) print(i_clust)
      for (j in 1:i-1) {
        j_clust = clusters[j]
        # if (i == 1) print(j_clust)
        new_df = rbind(new_df, df[which(df$df1_cluster == i_clust & df$df2_cluster == j_clust),])
      }
    }
    df = new_df
  }
  
  # Color for text label in heatmap
  df$id = rownames(df)
  df$ovlp_col = df$ovlp > mean(df$ovlp)
  df$pct_col = df$pct > mean(df$pct)
  df$ovlp_same_dir_col = df$ovlp_same_dir > mean(df$ovlp_same_dir)
  df$pct_same_dir_col = df$pct_same_dir > mean(df$pct_same_dir)
  
  # Find Max's
  df$ovlp_best = df$ovlp
  df$pct_best  = df$pct
  df$ovlp_same_dir_best = df$ovlp_same_dir
  df$pct_same_dir_best  = df$pct_same_dir
  for (cluster in unique(df$df1_cluster)) {
    # Find Ovlp Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$ovlp)]
    df$ovlp_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Pct Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$pct)]
    df$pct_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Ovlp Same Direction Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$ovlp_same_dir)]
    df$ovlp_same_dir_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Pct Same Direction Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$pct_same_dir)]
    df$pct_same_dir_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
  }
  
  
  if (any(sign(all_logFC) == -1)) {
    print("Up and downregulated DEGs analyzed")
    png1_name = paste(filepath, filename, "_ovlp_same_dir.png", sep="")
    png2_name = paste(filepath, filename, "_best_guess_same_dir.png", sep="")
    png3_name = paste(filepath, filename, "_pct_same_dir.png", sep="")
    png4_name = paste(filepath, filename, "_pct_best_guess_same_dir.png", sep="")
    png5_name = paste(filepath, filename, "_mag_dend_same_dir.png", sep="")
    png6_name = paste(filepath, filename, "_dend_same_dir.png", sep="")
    pdf3_name = paste(filepath, filename, "_pct_same_dir.pdf", sep="")
    
    png1_title = paste("DEGs in Common w/ Same Sign b/w Clusters")
    png2_title = paste("Best Guess of DEGs w/ Same Sign")
    png3_title = paste("% DEGs w/ Same Sign in Common Clusters", with_correction)
    png4_title = paste("% Best Guess of DEGs w/ Same Sign", with_correction)
    
    df$ovlp = df$ovlp_same_dir
    df$ovlp_col = df$ovlp_same_dir_col
    df$ovlp_best = df$ovlp_same_dir_best
    df$pct = df$pct_same_dir
    df$pct_col = df$pct_same_dir_col
    df$pct_best = df$pct_same_dir_best
  } else {
    print("Only Upregulated Genes")
    png1_name = paste(filepath, filename, "_ovlp.png", sep="")
    png2_name = paste(filepath, filename, "_best_guess.png", sep="")
    png3_name = paste(filepath, filename, "_pct.png", sep="")
    png4_name = paste(filepath, filename, "_pct_best_guess.png", sep="")
    png5_name = paste(filepath, filename, "_mag_dend.png", sep="")
    png6_name = paste(filepath, filename, "_dend.png", sep="")
    
    png1_title = paste("DEGs in Common b/w Clusters")
    png2_title = paste("Best Guess")
    png3_title = paste("% DEGs in Common b/w Clusters", with_correction)
    png4_title = paste("% Best Guess", with_correction)
  }
  
  df$df1_cluster = factor(df$df1_cluster, levels = rev(unique(df$df1_cluster)))
  df$df2_cluster = factor(df$df2_cluster, levels = unique(df$df2_cluster))
  print("Right before plotting")
  
  # Plot 1 - Ovlp
  png(png1_name, width = 250*length(dfs)+50, height = 250*length(dfs), unit = "px", res = 110)
  p = ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp)) + geom_tile() + scale_fill_viridis(discrete=FALSE) +  ggtitle(png1_title) + guides(color = FALSE) + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  if (labels)
    p = p + geom_text(aes(label=ovlp, color=ovlp_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000"))
  if (! xlab)
    p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  print(p)
  dev.off()
  print("finished plot 1")
  
  # Plot 2 - Ovlp Best Guess
  png(png2_name, width = 250*length(dfs)+50, height = 250*length(dfs), unit = "px", res = 110)
  p = ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + ggtitle(png2_title) + guides(color = FALSE) + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  if (labels)
    p = p + geom_text(data=subset(df, ovlp_same_dir_best > 0), aes(label=ovlp_same_dir_best, color=ovlp_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000"))
  if (! xlab)
    p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  print(p)
  dev.off()
  print("finished plot 2")
  
  # Plot 3 - Pct
  png(png3_name,  width = 250*length(dfs)+50, height = 250*length(dfs), unit = "px", res = 110)
  p = ggplot(df, aes(df1_cluster, df2_cluster, fill=pct)) + geom_raster() + scale_fill_viridis(discrete=FALSE) + ggtitle(png3_title) + guides(color = FALSE) + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  if (labels)
    p = p + geom_text(aes(label=format(round(pct, 1), nsmall = 1), color=pct_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) 
  if (! xlab)
    p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  print(p)
  dev.off()
  pdf(pdf3_name,  width = 2.5*length(dfs)+.50, height = 2.5*length(dfs), version = "1.6", bg = "white")
  print(p)
  dev.off()
  print("finished plot 3")
  
  # Plot 4 - Pct Best Guess
  png(png4_name,  width = 250*length(dfs)+50, height = 250*length(dfs), unit = "px", res = 110)
  p = ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_same_dir_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + ggtitle(png4_title) + guides(color = FALSE) + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  if (labels)
    p = p + geom_text(data=subset(df, pct_same_dir_best > 0), aes(label=format(round(pct_same_dir_best, 1), nsmall = 1), color=pct_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000"))
  if (! xlab)
    p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  print(p)
  dev.off()
  print("finished plot 4")
  
  if (dendrogram) {
    # Prepare the Data
    all_genes = all_genes[which(duplicated(all_genes))] # only keep those which are found in comparisons twice. Reduces number of genes
    all_genes = unique(all_genes)
    dend_mat = matrix(, nrow=length(all_genes), ncol=length(all_clusters), dimnames = list(all_genes, all_clusters))
    print(paste("Creating Dendrogram Matrix of size", nrow(dend_mat), "x", ncol(dend_mat)))
    for (i in 1:length(dfs)) {
      for (i_clust in unique(dfs[[i]]$cluster)) {
        i_clust_df = dfs[[i]][which(dfs[[i]]$cluster == i_clust),]
        i_clust_df = i_clust_df[!duplicated(i_clust_df$gene),]
        gene_in_all_genes = i_clust_df$gene[which(i_clust_df$gene %in% all_genes)]
        clust_name = paste(samples[[i]], i_clust)
        dend_mat[gene_in_all_genes, clust_name] = i_clust_df$avg_logFC[which(i_clust_df$gene %in% all_genes)]
      }
    }
    
    na_ind = which(is.na(dend_mat))
    dend_mat[na_ind] = 0
    # dend_mat=dend_mat[1:2000,]
    
    Colors=brewer.pal(11,"Spectral")
    # my_breaks = c(seq(-10, -2, length=100), seq(-2,-0.1,length=100), seq(-0.1, 0.1, length=100), seq(0.1,2,length=100), seq(2,10,length=100))
    my_palette1=colorRampPalette(Colors)(n = 299)
    my_palette <- colorRampPalette(c("#ff4b5c", "#FFFFFF", "#056674"))(n = 299)
    print("Plotting the magnitude dendrogram")
    par(mar=c(10, 4.1, 4.1, 2.1))
    png(png5_name, width = 300*length(dfs)+50, height = 300*length(dfs), unit = "px", res = 120)
    heatmap.2(dend_mat, scale = "none", dendrogram = "both", col = my_palette1, trace = "none", margins=c(10,5), srtCol=45, symm=F,symkey=F,symbreaks=T, main="Magnitude of Expression")
    dev.off()
    print("finished magnitude dendrogram")
    
    dend_mat = sign(dend_mat)
    dend_mat[na_ind] = 0
    print("Plotting the direction dendrogram")
    par(mar=c(10, 4.1, 4.1, 2.1))
    png(png6_name, width = 300*length(dfs)+50, height = 300*length(dfs), unit = "px", res = 120)
    heatmap.2(dend_mat, scale = "none", dendrogram = "both", col = my_palette1, trace = "none", margins=c(10,5), srtCol=45, symm=F,symkey=F,symbreaks=T, main="Direction of Expression")
    dev.off()
    print("finished direction dendrogram")
  }
  
  
  return(list(df, dend_mat))
}

heatmapComparison <- function(df1, df2, df1_sample, df2_sample, filename, filepath) {
  # Input: 2 dataframes that are output of Seurat FindAllMarkers
  #        The samples or whatever you want to name those two dataframes
  #        Base file name for png output
  #        Filepath for png output
  # Output: 4 png files. A comparison of DEGs in the two input dataframes.
  df1_clusters = unique(as.vector(df1$cluster))
  df2_clusters = unique(as.vector(df2$cluster))
  df1_num_clusters <- length(df1_clusters)
  df2_num_clusters <- length(df2_clusters)
  df <- data.frame()
  gene_df = data.frame()
  for (i in 1:df1_num_clusters) {
    for (j in 1:df2_num_clusters) {
      df1_cluster <- df1[which(df1$cluster == df1_clusters[i]),]
      df2_cluster <- df2[which(df2$cluster == df2_clusters[j]),]
      df1_cluster = df1_cluster[!duplicated(df1_cluster$gene),]
      df2_cluster = df2_cluster[!duplicated(df2_cluster$gene),]
      
      ovlp_genes = unique(df2_cluster$gene[which(df2_cluster$gene %in% df1_cluster$gene)])
      ovlp_genes = ovlp_genes[which(! is.na(ovlp_genes))]
      ovlp = length(unique(ovlp_genes))
      df1_sign = sign(df1_cluster$avg_logFC[which(df1_cluster$gene %in% ovlp_genes)])
      df2_sign = sign(df2_cluster$avg_logFC[which(df2_cluster$gene %in% ovlp_genes)])
      names(df1_sign) = df1_cluster$gene[which(df1_cluster$gene %in% ovlp_genes)]
      names(df2_sign) = df2_cluster$gene[which(df2_cluster$gene %in% ovlp_genes)]
      if (length(df1_sign) != length(df2_sign)) { print("Error in signs."); return(); }
      df2_sign = df2_sign[names(df1_sign)]
      ovlp_same_dir_genes = unique(names(df1_sign[which(df1_sign == df2_sign)]))
      ovlp_same_dir = length(ovlp_same_dir_genes)
      
      total_ovlp = 2*ovlp
      total_ovlp_same_dir = 2*ovlp_same_dir
      if ("correction_factor" %in% colnames(df2_cluster) && "correction_factor" %in% colnames(df1_cluster)) {
        total_ovlp = ovlp*df1_cluster$correction_factor[1] + ovlp*df2_cluster$correction_factor[1]
        total_ovlp_same_dir = ovlp_same_dir*df1_cluster$correction_factor[1] + ovlp_same_dir*df2_cluster$correction_factor[1]
        with_correction = "w/ Correction for Gene Conversion"
        # print(paste0("Correction factor found in ", df1_sample, "_", df1_clusters[i], " and ", df2_sample, "_", df2_clusters[j]))
        # print(colnames(df1_cluster))
        # print(colnames(df2_cluster))
      }
      
      pct = (total_ovlp / (nrow(df1_cluster) + nrow(df2_cluster))) * 100
      pct_same_dir = (total_ovlp_same_dir / (nrow(df1_cluster) + nrow(df2_cluster))) * 100
      
      # Check if pct is greater than 100
      if (pct_same_dir > 100) {
        print("Error pct ovlp > 100")
        print(paste0(df1_sample, "_", df1_clusters[i]))
        print(paste0(df2_sample, "_", df2_clusters[j]))
        print(ovlp_same_dir)
        print(nrow(i_clust_df))
        print(nrow(j_clust_df))
        print(total_ovlp_same_dir)
      }
      
      df <- rbind(df, t(c(df1_clusters[i], df2_clusters[j], ovlp, pct, ovlp_same_dir, pct_same_dir)))
      
      new_gene_df_rows = data.frame(rep(paste(df1_clusters[i]), ovlp_same_dir), rep(paste(df2_clusters[j]), ovlp_same_dir), sort(ovlp_same_dir_genes), df1_sign[sort(ovlp_same_dir_genes)])
      colnames(new_gene_df_rows) = c(paste(df1_sample, "Cluster"), paste(df2_sample, "Cluster"), "Genes In Common", "Direction of Regulation")
      gene_df = rbind(gene_df, new_gene_df_rows)
    }
  }
  
  colnames(df) <- c("df1_cluster", "df2_cluster", "ovlp", "pct", "ovlp_same_dir", "pct_same_dir")
  df$df1_cluster = factor(df$df1_cluster, levels = df1_clusters)
  df$df2_cluster = factor(df$df2_cluster, levels = df2_clusters)
  df$ovlp = as.numeric(as.vector(df$ovlp))
  df$pct = as.numeric(as.vector(df$pct))
  df$ovlp_same_dir = as.numeric(as.vector(df$ovlp_same_dir))
  df$pct_same_dir = as.numeric(as.vector(df$pct_same_dir))
  
  # Sort Clusters Numerically if Possible
  # if (! any(is.na(as.numeric(df$df1_cluster))))
  #   df$df1_cluster = as.numeric(df$df1_cluster)
  # if (! any(is.na(as.numeric(df$df2_cluster))))
  #   df$df1_cluster =as.numeric(df$df2_cluster)
  
  
  # Color for text label in heatmap
  df$id = rownames(df)
  df$ovlp_col = df$ovlp > mean(df$ovlp)
  df$pct_col = df$pct > mean(df$pct)
  df$ovlp_same_dir_col = df$ovlp_same_dir > mean(df$ovlp_same_dir)
  df$pct_same_dir_col = df$pct_same_dir > mean(df$pct_same_dir)
  
  # Find Max's
  df$ovlp_best = df$ovlp
  df$pct_best  = df$pct
  df$ovlp_same_dir_best = df$ovlp_same_dir
  df$pct_same_dir_best  = df$pct_same_dir
  for (cluster in unique(df$df1_cluster)) {
    # Find Ovlp Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$ovlp)]
    df$ovlp_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Pct Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$pct)]
    df$pct_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Ovlp Same Direction Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$ovlp_same_dir)]
    df$ovlp_same_dir_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Pct Same Direction Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$pct_same_dir)]
    df$pct_same_dir_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
  }
  
  # Plot 1 - Ovlp
  if (any(sign(df1$avg_logFC) == -1) || any(sign(df2$avg_logFC) == -1)) {
    png(paste(filepath, filename, "_ovlp_same_dir.png", sep=""), width = df1_num_clusters*100, height = df2_num_clusters*100, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp_same_dir)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=ovlp_same_dir, color=ovlp_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("DEGs in Common w/ Same Sign b/w", df1_sample, "&", df2_sample,  "Clusters")) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()) + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
    dev.off()
  } else {
    png(paste(filepath, filename, "_ovlp.png", sep=""),  width = df1_num_clusters*100, height = df2_num_clusters*100, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=ovlp, color=ovlp_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("DEGs in Common b/w", df1_sample, "&", df2_sample,  "Clusters")) + guides(color = FALSE) + theme_classic() + theme(line = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
    dev.off()
  }
  print("finished 1")

  # Plot 2 - Ovlp Best Guess
  if (any(sign(df1$avg_logFC) == -1) || any(sign(df2$avg_logFC) == -1)) {
    png(paste(filepath, filename, "_best_guess_same_dir.png", sep=""), width = df1_num_clusters*100, height = df2_num_clusters*100, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp_same_dir_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(data=subset(df, ovlp_same_dir_best > 0), aes(label=ovlp_same_dir_best, color=I("#000000"))) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("Best Guess of DEGs w/ Same Sign b/w", df1_sample, "&", df2_sample)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
    dev.off()
  } else {
    png(paste(filepath, filename, "_best_guess.png", sep=""), width = df1_num_clusters*100, height = df2_num_clusters*100, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(data=subset(df, ovlp_best > 0), aes(label=ovlp_best, color=I("#000000"))) +  xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("Best Guess b/w", df1_sample, "&", df2_sample)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
    dev.off()
  }
  print("finished 2")

  # Plot 3 - Pct
  if (any(sign(df1$avg_logFC) == -1) || any(sign(df2$avg_logFC) == -1)) {
    png(paste(filepath, filename, "_pct_same_dir.png", sep=""), width = df1_num_clusters*100, height = df2_num_clusters*100, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_same_dir)) + geom_raster() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=format(round(pct_same_dir, 1), nsmall = 1), color=pct_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% DEGs w/ Same Sign in Common b/w", df1_sample, "&", df2_sample,  "Clusters", with_correction)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
    dev.off()
    pdf(paste(filepath, filename, "_pct_same_dir.pdf", sep=""), width = df1_num_clusters*1.00 + 0.5, height = df2_num_clusters*1.00, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_same_dir)) + geom_raster() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=format(round(pct_same_dir, 1), nsmall = 1), color=pct_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% DEGs w/ Same Sign in Common b/w", df1_sample, "&", df2_sample,  "Clusters", with_correction)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
    dev.off()
  } else {
    png(paste(filepath, filename, "_pct.png", sep=""), width = df1_num_clusters*100, height = df2_num_clusters*100, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=format(round(pct, 1), nsmall = 1), color=pct_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% DEGs in Common b/w", df1_sample, "&", df2_sample,  "Clusters", with_correction)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
    dev.off()
  }
  print("finished plot 3")
  
  # Plot 4 - Pct Best Guess
  if (any(sign(df1$avg_logFC) == -1) || any(sign(df2$avg_logFC) == -1)) {
    png(paste(filepath, filename, "_pct_best_guess_same_dir.png", sep=""), width = df1_num_clusters*100, height = df2_num_clusters*100, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_same_dir_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(data=subset(df, pct_same_dir_best > 0), aes(label=format(round(pct_same_dir_best, 1), nsmall = 1), color=I("#000000"))) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% Best Guess of DEGs w/ Same Sign b/w", df1_sample, "&", df2_sample)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
    dev.off()
  } else {
    png(paste(filepath, filename, "_pct_best_guess.png", sep=""), width = df1_num_clusters*100, height = df2_num_clusters*100, unit = "px", res = 100)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(data=subset(df, pct_best > 0), aes(label=format(round(pct_best, 1), nsmall = 1), color=I("#000000"))) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% Best Guess b/w", df1_sample, "&", df2_sample)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
    dev.off()
  }
  print("finished plot 4")
  
  return(list(df, gene_df))
}

downsampleObj = function(obj, run = 1) {
  #' Downsample the counts of the Seurat Object to the cell with the lowest number of counts.
  #' 
  #' @param obj Seurat object to downsample
  #' @param run the run number, also used as the seed for random sampling
  #' @return the downsampled object
  
  new_mat = downsample(obj, rownames(obj), run)
  obj = SetAssayData(obj, slot = "counts", new.data = new_mat)
  obj = NormalizeData(obj)
  obj = ScaleData(obj)
  
  return(obj)
}

downsample <- function(combined, marker_genes, run) {
  set.seed(run)
  min_trans <- min(combined$nCount_RNA)
  gene_names <- rownames(combined@assays$RNA@counts)
  # new_matrix <- matrix(, nrow = nrow(combined@assays$RNA@counts), ncol = ncol(combined@assays$RNA@counts), dimnames = list(gene_names, colnames(combined@assays$RNA@counts)))
  # new_new_matrix <- matrix(, nrow=nrow(combined@assays$RNA@counts))
  marker_matrix  <- matrix(, nrow=length(marker_genes), ncol = ncol(combined@assays$RNA@counts), dimnames = list(marker_genes, colnames(combined@assays$RNA@counts)))
  i <- 0
  for (cell in colnames(combined@assays$RNA@counts)) {
    # i <- i + 1
    # if (i%%500 == 1) {
    #   print(cell)
    # }
    # start.time <- Sys.time()
    
    trans_names <- rep(gene_names, combined@assays$RNA@counts[,cell])
    ran_trans_names <- sample(trans_names, min_trans)
    ran_trans_names <- ran_trans_names[which(ran_trans_names %in% marker_genes)]
    ran_df <- as.data.frame(table(ran_trans_names))
    zero_gene_names <- marker_genes[which(! marker_genes %in% ran_trans_names)]
    zero_df <- setNames(data.frame(zero_gene_names <- zero_gene_names, Freq <- rep(0, length(zero_gene_names))), c("ran_trans_names", "Freq"))
    ran_df <- rbind(ran_df, zero_df)
    rownames(ran_df) <- ran_df$ran_trans_names
    ran_df <- ran_df[marker_genes,2]
    # new_matrix[,cell] <-as.matrix(ran_df)
    
    
    marker_matrix[,cell] <- as.matrix(ran_df)
    # new_new_matrix <- cbind(new_new_matrix, as.matrix(ran_df))
    
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # print(time.taken)
  }
  
  return(marker_matrix)
}

pickNewCells <- function(combined, num_clusters, num_cells) {
  new_cells <- c()
  for (i in 1:num_cells) {
    ran_cluster <- sample(0:num_clusters, 1)
    this_cells <- names(combined$seurat_clusters[which(combined$seurat_clusters == ran_cluster)])
    new_cells <- c(new_cells, sample(this_cells,1))
  }
  
  return(new_cells)
}

shuffleClusters <- function(combined) {
  # The selection process for a new cluster should be as follows:
  # 1. Pick a random cluster 0-40
  # 2. Pick a random cell from that cluster to be a part of the new cluster
  # This means that the new data set would likely have duplicate cells
  new_cells <- lapply(0:num_clusters, function(x) c())
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  for (i in 0:num_clusters) {
    num_cells <- length(combined$seurat_clusters[which(combined$seurat_clusters == i)])
    new_cells[[i+1]] <- pickNewCells(combined, num_clusters, num_cells)
  }
  
  return(new_cells)
}