
# Load libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(cowplot)
library(GO.db)

setwd('/Volumes/ftp/Paula Sole/mixOmics/AR/GOEA')

# Load data
data <- read.csv("/Volumes/ftp/Paula Sole/mixOmics/AR/GOEA/network_AR_0.csv", header=TRUE, sep=",")
Genes <- toupper(data$label)
Genes <- as.character(Genes)


# Load GO annotations
GOTERM <- read.table("~/Desktop/GOTERM.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(GOTERM) <- c("GO", "GeneID")
GOTERM$GeneID <- toupper(GOTERM$GeneID)  # Ensure consistency in case sensitivity

# Replace GO IDs with a list of GO descriptions
GO_ids <- unique(GOTERM$GO)
GO_descriptions <- AnnotationDbi::select(GO.db, 
                                         keys = GO_ids, 
                                         keytype = "GOID", 
                                         columns = c("GOID", "TERM"))
head(GO_descriptions)

# Ensure no NA values are in TERM
GO_descriptions <- GO_descriptions[!is.na(GO_descriptions$TERM), ]

# Extract all genes (background)
background_genes <- unique(GOTERM$GeneID)

#Perform GO enrichment analysis
gse <- enricher(gene = Genes,
                pvalueCutoff = 0.05,
                pAdjustMethod = "none",
                universe = background_genes,
                minGSSize = 30,
                maxGSSize = 500,
                TERM2GENE = GOTERM
                #TERM2NAME = GO_descriptions
                )
  
# Check if enrichment analysis returned results
if (is.null(gse) || nrow(gse@result) == 0) {
  stop("No significant enrichment results.")
}
  
# Save plots to PDF
pdf(file = "GOEnrichment_plots_AR_GO.pdf")

# Barplot
if (nrow(gse@result) > 0) {
  print(barplot(gse,
                x = "GeneRatio",
                showCategory = 20,
                font.size = 9) + 
          ggtitle("Barplot of GO Enrichment"))
  
  # Dotplot
  print(dotplot(gse, 
                showCategory = 20,
                font.size = 8) + 
          ggtitle("Dotplot of GO Enrichment"))
}

dev.off()


# Save cnetplot to PDF
pdf(file = "Cnetplot_GOEnrichment_AR_descriptions.pdf")

cnetplot(gse,
         layout = igraph::layout_with_kk,
         showCategory = 5,
         color_category = "#E5C494",
         size_category = 1,
         color_item = "#B3B3B3",
         size_item = 1,
         color_edge = "grey",
         size_edge = 0.5,
         node_label = "all",
         foldChange = NULL,
         hilight = "none",
         hilight_alpha = 0.3)

dev.off()
