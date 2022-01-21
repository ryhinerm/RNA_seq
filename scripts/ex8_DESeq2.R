library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(BiocParallel)
library(org.Rn.eg.db) # Rat's Annotationdbi
register(MulticoreParam(4)) # Change this based on your computer core count

setwd("<Insert here your working directory>")

# Setting up color profiles from colorbrewer
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_1 <- "Neuropil_Poly"
sample_2 <- "Somata_Poly"
no_of_reps <- 3

sample_column <- c(rep(sample_1, no_of_reps),
                   rep(sample_2, no_of_reps))

run_column <- c(paste(sample_1, "1", sep = "_"),
                paste(sample_1, "2", sep = "_"),
                paste(sample_1, "3", sep = "_"),
                paste(sample_2, "1", sep = "_"),
                paste(sample_2, "2", sep = "_"),
                paste(sample_2, "3", sep = "_"))

rep_column <- c("A", "B", "C",
                "A", "B", "C")

samples_df <- data.frame(sample_column,
                         run_column,
                         rep_column)

colnames(samples_df) <- c("sample", "run", "rep")

samples_df$condition <- factor(rep(c(sample_1, sample_2), each = no_of_reps))

rownames(samples_df) <- samples_df$run

# Load count data
featurecount_data <- read.table("CDS_counts_processed.txt", header = TRUE, row.names = 1)

# Change colnames
# Make sure the column order in featurecount_data matches samples_df !!!
colnames(featurecount_data) <- rownames(samples_df)

#Import as DESeqDataSet (dds)
dds <- DESeqDataSetFromMatrix(countData = featurecount_data,
                              colData = samples_df,
                              design = ~ condition)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Factor levels
# WT columns should be first followed by KO/treatment
dds$condition <- factor(c(rep(sample_1, no_of_reps),
                          rep(sample_2, no_of_reps)),
                        levels = c(sample_1,
                                   sample_2))

# Differential expression analysis
dds <- DESeq(dds)

# colData(dds) # to check whether names are correct

################################################################################
################################################################################
# QC
################################################################################
################################################################################

# Log transformation for data quality assessment
rld <- rlog(dds, blind = FALSE)

# Sample distance matrix
sampleDists <- as.matrix(dist(t(assay(rld))))
pdf("QC_sample_distance_matrix_CDS.pdf")
heatmap.2(as.matrix(sampleDists),
          key = T,
          trace = "none",
          col = colorpanel(100, "#2b8cbe", "#7bccc4", "white"),
          ColSideColors = mycols[dds$condition],
          RowSideColors = mycols[dds$condition],
          margin = c(10, 10), main = "Sample Distance Matrix")
dev.off()

# Count matrix heatmap
select <- order(rowMeans(counts(dds,normalized = TRUE)))
df <- as.data.frame(colData(dds)[ , c("condition","rep")])

pdf("QC_count_matrix_CDS.pdf")
pheatmap(assay(rld)[select,],
         cluster_rows = FALSE,
         show_rownames=FALSE,
         cluster_cols = FALSE,
         annotation_col = df)
dev.off()

# PCA plot

pdf("QC_PCA_CDS.pdf")
pcaData <- plotPCA(rld, intgroup = c("condition", "rep"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = rep)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()
################################################################################
################################################################################

################################################################################
################################################################################

# Using contrast so that the following code can be scaled with increasing number of samples
res <- results(dds, contrast = c("condition", sample_1, sample_2), alpha = 0.05)

# Adding gene names using org.Rn.eg.db
# Source: http://bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html
# Also: https://support.bioconductor.org/p/66288/
# This function takes a list of IDs as first argument and their key type as the second argument.
# The third argument is the key type we want to convert to, the fourth is the AnnotationDb object to use.
# Finally, the last argument specifies what to do if one source ID maps to several target IDs:
# should the function return an NA or simply the first of the multiple IDs

convertIDs <- function( ids, from, to, db, ifMultiple = c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

# # Check columns in the database that you want to add:
columns(org.Rn.eg.db)

# Actual adding of the column
res$GeneID <- row.names(res)
res$gene_symbol <- convertIDs(row.names(res), "ENSEMBL", "SYMBOL", org.Rn.eg.db)

summary(res)

res_df <- as.data.frame(res)

# Which genes are translationally deregulated:
res_df$regulation_level <- ifelse((res_df$log2FoldChange > 0.5 & res_df$padj < 0.05), "Upregulated",
                                  ifelse((res_df$log2FoldChange < - 0.5 & res_df$padj < 0.05),
                                         "Downregulated", "Unchanged"))

write.table(res_df,
            file = "DESeq2_res.csv",
            sep = ",",
            row.names = F,
            col.names = T,
            quote = F)

res_df$regulation_level <- factor(res_df$regulation_level, levels = c("Upregulated", "Downregulated", "Unchanged"))

res_df <- res_df[!is.na(res_df$padj), ]

# Plot

pdf("Volcano_plot.pdf", width = 4, height = 5)
ggplot(res_df,
       aes(x = -log10(padj),
           y = log2FoldChange,
           color = regulation_level)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#404788FF", "#73D055FF", "#999999")) +
  xlab("-log10(adjusted p-value)") +
  ylab("Log2 fold change") +
  labs(color = "Regulation level") +
  theme_bw()
dev.off()