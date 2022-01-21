library(grid)
library(gridExtra)
library(pathview)
library(AnnotationDbi)
library(clusterProfiler)
library(topGO)
library(ggplot2)

setwd("<Insert here your working directory>")

sample_name = "Neuropil_Poly_vs_Somata_Poly_only_padj"

df <- read.csv("DESeq2_res.csv", sep = ",", header = T)

rownames(df) <- df$GeneID
df <- df[order(df$padj), ]

# Define upregulated and downregulated genes based on padj value

genes_up <- which(df$padj < 0.05 & df$log2FoldChange > 0)
genes_down <- which(df$padj < 0.05 & df$log2FoldChange < 0)

all_genes_names <- rownames(df)

genes_up <- rownames(df)[genes_up]
genes_down <- rownames(df)[genes_down]

genelist_up <- factor(as.integer(all_genes_names %in% genes_up))
names(genelist_up) <- all_genes_names

genelist_down <- factor(as.integer(all_genes_names %in% genes_down))
names(genelist_down) <- all_genes_names

allGO2genes <- annFUN.org(whichOnto = "ALL",
                          feasibleGenes = NULL,
                          mapping = "org.Rn.eg.db",
                          ID = "ensembl")

GOdata_up_bp <- new("topGOdata",
                    ontology = "BP",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes,
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_up_mf <- new("topGOdata",
                    ontology = "MF",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes,
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_up_cc <- new("topGOdata",
                    ontology = "CC",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes, 
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_down_bp <- new("topGOdata",
                      ontology = "BP",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

GOdata_down_mf <- new("topGOdata",
                      ontology = "MF",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

GOdata_down_cc <- new("topGOdata",
                      ontology = "CC",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)


resultFis_up_bp <- runTest(GOdata_up_bp, statistic = "fisher")
resultFis_up_mf <- runTest(GOdata_up_mf, statistic = "fisher")
resultFis_up_cc <- runTest(GOdata_up_cc, statistic = "fisher")
resultFis_down_bp <- runTest(GOdata_down_bp, statistic = "fisher")
resultFis_down_mf <- runTest(GOdata_down_mf, statistic = "fisher")
resultFis_down_cc <- runTest(GOdata_down_cc, statistic = "fisher")

parse_tables <- function(GO_data, statistics)
{
  goEnrichment <- GenTable(GO_data, weightFisher = statistics, topNodes = 20)
  sub("< ", "", goEnrichment$weightFisher)
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$weightFisher <- as.numeric(sub("< ", "", goEnrichment$weightFisher))  
  goEnrichment
}

GOres_up_bp <- parse_tables(GOdata_up_bp, resultFis_up_bp)
GOres_up_mf <- parse_tables(GOdata_up_mf, resultFis_up_mf)
GOres_up_cc <- parse_tables(GOdata_up_cc, resultFis_up_cc)

GOres_down_bp <- parse_tables(GOdata_down_bp, resultFis_down_bp)
GOres_down_mf <- parse_tables(GOdata_down_mf, resultFis_down_mf)
GOres_down_cc <- parse_tables(GOdata_down_cc, resultFis_down_cc)


plot_GO <- function(GO_data, Ontology, Regulation, use_color) {
  GO_data$log_weightFisher <- (- log10(as.numeric(GO_data$weightFisher)))
  ggplot(GO_data, 
         aes(x = GO_data$log_weightFisher,
             y = GO_data$Term)) +
    geom_segment(aes(x = 0,
                     xend = GO_data$log_weightFisher,
                     y = GO_data$Term,
                     yend = GO_data$Term),
                 colour = use_color)  +
    geom_point(aes(size = GO_data$Significant),
               colour = use_color) +
    scale_size_area(name = "Gene counts") +
    xlab("Enrichment (- log10 Pvalue)") +
    ylab(Ontology) +
    ggtitle(Regulation) +
    scale_x_continuous() +
    theme_bw() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"))
}

plot_up_BP <- plot_GO(GOres_up_bp, "Biological Proccess", "TopGO Up (fisher's exact test)", "#404788FF")
plot_up_MF <- plot_GO(GOres_up_mf, "Molecular Function", "TopGO Up (fisher's exact test)", "#404788FF")
plot_up_CC <- plot_GO(GOres_up_cc, "Cellular Component", "TopGO Up (fisher's exact test)", "#404788FF")

# grid.arrange(grobs = list(p1,p2,p3))

plot_down_BP <- plot_GO(GOres_down_bp, "Biological Proccess", "TopGO Down (fisher's exact test)", "#73D055FF")
plot_down_MF <- plot_GO(GOres_down_mf, "Molecular Function", "TopGO Down (fisher's exact test)", "#73D055FF")
plot_down_CC <- plot_GO(GOres_down_cc, "Cellular Component", "TopGO Down (fisher's exact test)", "#73D055FF")

# grid.arrange(grobs = list(p1,p2,p3))

pdf(paste(sample_name, "Biological_Proccess_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_BP
dev.off()

pdf(paste(sample_name, "Molecular_Function_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_MF
dev.off()

pdf(paste(sample_name, "Cellular_Component_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_CC
dev.off()

pdf(paste(sample_name, "Biological_Proccess_TopGO_down_fisher.pdf", sep = "_"))
plot_down_BP
dev.off()

pdf(paste(sample_name, "Molecular_Function_TopGO_down_fisher.pdf", sep = "_"))
plot_down_MF
dev.off()

pdf(paste(sample_name, "Cellular_Component_TopGO_down_fisher.pdf", sep = "_"))
plot_down_CC
dev.off()