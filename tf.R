rna_path <- "C:/Users/miles/Downloads/brain/"
tf_df <- read.table("C:/Users/miles/Downloads/brain/brain_scripts/tf_shiny/data/Homo-sapiens_all.txt", stringsAsFactors = FALSE, header = FALSE)

all_tf <- unique(tf_df[,1])
test <- data.frame()
for (tf in all_tf) {
  this_rows <- tf_df[which(tf_df[,1] == tf),]
  test <- rbind(test, t(c(tf, nrow(this_rows))))
}
colnames(test) <- c("tf", "num_genes")
test$num_genes <- as.numeric(as.vector(test$num_genes))
ggplot(test, aes(x = "", y = num_genes)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(width = 0.38), alpha = 0.2)
