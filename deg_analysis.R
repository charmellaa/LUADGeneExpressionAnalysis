rm(list=ls())
options(stringAsFactors = F)

library(stringr)
library(pheatmap)


# 1, PREPARING THE DATA
dirRes <- "Results/"

if (!dir.exists(dirRes)) {
  dir.create(dirRes)
} else {
  print(paste("The directory", dirRes, "already exists"))
}

tumor <- "LUAD"

dirTumor <- paste0(dirRes, tumor, "/")
if (!dir.exists(dirTumor)) {
  dir.create(dirTumor)
} else {
  print(paste("The directory", dirTumor, "already exists"))
}

path <- "./TCGA_datasets/matrix"
filename_in <- paste0(path, "/matrix_RNAseq_", tumor, ".txt")

filename_DEG <- paste0(dirTumor, "DEG.txt")
filename_matrix_DEG <- paste0(dirTumor, "matrix_DEG.txt")
filename_heatmap <- paste0(dirTumor, "heatmap.pdf")

filename_list_normal <- paste0(dirTumor, "normal.txt")
filename_list_tumor <- paste0(dirTumor, "tumor.txt")


# parameters

percentile <- 0.1
threshold_fc <- 2
threshold_pval <- 0.05
paired <- TRUE

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, row.names = 1, quote = "", nrows = 10)
classes <- sapply(tmp, class)
tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, row.names = 1, quote = "", colClasses = classes)

genes <- rownames(tmp) # vector
patients <- colnames(tmp)

normalPatients <- grep('TCGA-\\w+-\\w+-1\\d', patients, value = TRUE)
cancerPatients <- grep('TCGA-\\w+-\\w+-0\\d', patients, value = TRUE)

normalPatients <- normalPatients[!duplicated(str_extract(normalPatients, "TCGA-\\w+-\\w+"))]
cancerPatients <- cancerPatients[!duplicated(str_extract(cancerPatients, "TCGA-\\w+-\\w+"))]

nameN <- str_extract(normalPatients, "TCGA-\\w+-\\w+")
nameC <- str_extract(cancerPatients, "TCGA-\\w+-\\w+")

common_patients <- intersect(nameN, nameC)
                            
normalPatients_ids <- unlist(lapply(common_patients, function(x) {grep(x, normalPatients, value = TRUE)}))
cancerPatients_ids <- unlist(lapply(common_patients, function(x) {grep(x, cancerPatients, value = TRUE)}))

data_normal <- tmp[,normalPatients_ids]
data_cancer <- tmp[,cancerPatients_ids]

data_all <- cbind(data_normal, data_cancer)
patients_ids <- colnames(data_all)

rm(tmp, patients, nameN, nameC, normalPatients, cancerPatients, common_patients, classes)



# 2. ANALYSIS

# delete obs with mean = 0
overall_mean <- apply(data_all, 1, mean)
i <- which(overall_mean == 0)
data_normal <- data_normal[-i,]
data_cancer <- data_cancer[-i,]
data_all <- data_all[-i,]
genes <- genes[-i]


# transform in logarithmic base 2
data_normal <- log2(data_normal+1)
data_cancer <- log2(data_cancer+1)
data_all <- log2(data_all+1)


# compute IQR
variation <- apply(data_all, 1, IQR, type = 5)

# filter
threshold_percentile <- quantile(variation, percentile)
i <- which(variation <= threshold_percentile)


data_normal <- data_normal[-i,]
data_cancer <- data_cancer[-i,]
data_all <- data_all[-i,]
genes <- genes[-i]

rm(i)


hist(variation,
     main = "IQR frequency distribution", breaks = 100,
     xlab = "IQR value", ylab = "Frequency", col = "cyan")

abline(v = threshold_percentile, lty = 2, lwd = 4, col = "pink")

# fold change
logFoldChange <- rowMeans(data_cancer) - rowMeans(data_normal)

hist(logFoldChange, main = "FC (log) frequency distribution", breaks = 100,
     xlab = "log Fold Change", ylab = "Frequency", col = "red")
abline(v = c(-log2(threshold_fc), log2(threshold_fc)), lty = 2, lwd = 4, col = "blue")


i <- which(abs(logFoldChange) < log2(threshold_fc))

if (length(i)>0) { 
  data_normal <- data_normal[-i,]
  data_cancer <- data_cancer[-i,]
  data_all <- data_all[-i,]
  genes <- genes[-i]
  logFoldChange <- logFoldChange[-i]
}

rm(i)

# p-value

colNormal <- ncol(data_normal)
colCancer <- ncol(data_cancer)

pvalue <- apply(data_all, 1, function(x) {
  t.test(x[1:colNormal], x[(colNormal+1):(colCancer+colNormal)], paired = paired)$p.value
})

pval_adjusted <- p.adjust(pvalue, method = "fdr")

hist(pval_adjusted, main = "Adjusted p-val Frequency distribution", breaks = 100,
     xlab = "adjusted p-value", ylab = "Frequency", col = "red")
abline(v = threshold_pval, lty = 2, lwd = 4, col = "blue")

i <- which(pval_adjusted > threshold_pval)

# hist(pval_adjusted)

# as before check ind length

if (length(i)>0) {
  data_normal <- data_normal[-i,]
  data_cancer <- data_cancer[-i,]
  data_all <- data_all[-i,]
  genes <- genes[-i]
  logFoldChange <- logFoldChange[-i]
  pvalue <- pvalue[-i]
  pval_adjusted <- pval_adjusted[-i]
}

rm(i)


# export results

direction <- ifelse(logFoldChange>0, "UP", "DOWN")

DEG <- data.frame(str_split_fixed(genes, "\\|", 2))
colnames(DEG) <- c("geneSymbol", "ensembl_id")

results <- data.frame(genes = DEG$geneSymbol, ensembl_id = DEG$ensembl_id, logFoldChange = logFoldChange, pvalue = pvalue, pvalue_adj = pval_adjusted, direction = direction)

# we can also order results based on an argument
# results <- results[order(results$logFC, decreasing = T),]

write.table(results, file = filename_DEG, row.names = F, sep = "\t", quote = F)
write.table(data_all, file = filename_matrix_DEG, row.names = T, col.names = NA, sep = "\t", quote = F)
write.table(normalPatients_ids, filename_list_normal, sep = "\t", col.names = F, row.names = F, quote = F)
write.table(cancerPatients_ids, filename_list_tumor, sep = "\t", col.names = F, row.names = F, quote = F)

# plots

plot(logFoldChange, -log10(pval_adjusted), 
     main = "Volcano plot", 
     # xlim and ylim depend on the datasets
     xlim = c(-7, 7),
     ylim = c(0,40),
     xlab = "log2 fold change", 
     ylab = "-log10, p-value")

abline(h = -log10(threshold_pval), 
       lty = 2, # line type, 2 = dashed
       lwd = 4,
       col = "red")

abline(v = c(-log2(threshold_fc), log2(threshold_fc)),
       lty = 2, # line type, 2 = dashed
       lwd = 4,
       col = "blue")

# boxplot on most upregulated gene
i <- which.max(results$logFoldChange) 
gene_id <- results$genes[i]

df = data.frame(normal = t(data_normal[i,]),
                cancer = t(data_cancer[i,]),
                row.names = NULL)

colnames(df) <- c(paste(gene_id, "normal"),
                  paste(gene_id, "cancer"))

boxplot(df,
        main = paste0(gene_id, ",", "adjusted p-value = ",
                      format(pval_adjusted[i], digits = 2)),
        notch = T, # boxplot with notch, or not if F 
        # notch is the interval of confidence around the median value
        # If the notches of two boxplots (representing different groups or categories) do not overlap, 
        #it suggests that there is strong evidence to suggest that the medians of the two groups are significantly different. 
        #On the other hand, if the notches overlap, it suggests that there is less evidence to suggest a difference in medians between the groups.
        ylab = "Gene expression value",
        col = c("green", "orange"),
        pars = list(boxwex = 0.3, staplewex = 0.6)
)

# boxplot on most downregulated gene
rm(i)
i <- which.min(results$logFoldChange) 
gene_id <- results$genes[i]

df = data.frame(normal = t(data_normal[i,]),
                cancer = t(data_cancer[i,]),
                row.names = NULL)

colnames(df) <- c(paste(gene_id, "normal"),
                  paste(gene_id, "cancer"))


boxplot(df,
        main = paste0(gene_id, ",", "adjusted p-value = ",
                      format(pval_adjusted[i], digits = 2)),
        notch = T,
        ylab = "Gene expression value",
        col = c("green", "orange"),
        pars = list(boxwex = 0.3, staplewex = 0.6)
)

rm(i)


# PIE CHART

count <- table(results$direction)
pie(count,
    labels = paste0(names(count), " ", round(100 * count/sum(count), 2), "%"),
    col = c("blue", "gold"))


# HEATMAP
# colored map
# clustering on the row and on the columns
# distance matrix for example correlation
# method of clustering (linkage)
test <- grepl('TCGA-\\w+-\\w+-1\\d', colnames(data_all))

condition <- ifelse(test, "normal", "cancer")

annotation <-  data.frame(condition = condition)
rownames(annotation) <- colnames(data_all)

vect_color <- c("green", "orange")
names(vect_color) <- unique(condition)

annotation_colors <- list(condition = vect_color)

pheatmap(data_all, scale = "row",
         border_color = NA,
         cluster_cols = T,
         cluster_rows = T,
         clustering_distance_rows = "correlation", #genes highly correlated close to each other
         clustering_distance_cols = "correlation",
         clustering_method = "centroid", #linkage, two sets of points how are they near?
         annotation_col = annotation, #normal green, cancer orange, for the column
         annotation_colors = annotation_colors,
         color = colorRampPalette(colors = c("blue3", "blue4", "brown4", "yellow2", "yellow"))(100), #from low expression genes blue to high yellow
         show_rownames = F,
         show_colnames = F,
         cutree_cols = 2, #first two main clustering cutted
         cutree_rows = 2,
         width = 10, height = 10,
         filename = filename_heatmap
)

