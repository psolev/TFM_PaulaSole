
# Title: DIABLO
# Created: November 2024
# Author: Paula Solé

# ---------------------------- Loading packages and data ------------------------------------
# Load packages
library(mixOmics)
library(gridExtra)
library(igraph)
library(readxl)
library(openxlsx)

# Clear your workplace
rm(list = ls())

# Change working directory
directory <- '//150.128.81.153/ftp/Paula Sole/chemical priming' # Change with your working directory
setwd(directory)

# Load data (samples in rows and features in columns)

# Metabolomics data
met_archive <- "Annotated_Met_QP.xlsx" # Change with your archive name 
met_sheet <- "Sheet1" 
metabolomics <-  read_excel(met_archive, sheet = met_sheet)
dim(metabolomics) # Check the dimension of the dat
metabolomics <- t(metabolomics) # Transpose the data
colnames(metabolomics) <- metabolomics[1,] # Set the metabolites as the column names
metabolomics <- metabolomics[-1,]
dim(metabolomics) # Check the dimension of the dat

# Transcriptomics data
trans_archive <- "processed_data.csv" # Change with your archive name 
transcriptomics <- read.csv(trans_archive, header=TRUE, sep=",")
dim(transcriptomics) # Check the dimension of the data
transcriptomics <- t(transcriptomics) # Transpose the data
colnames(transcriptomics) <- transcriptomics[1,] # Set the gene ID as the column names
transcriptomics <- transcriptomics[-(1:3),]
dim(transcriptomics) # Check the dimension of the data

# To change the gene names
#genes_annotation <- read.table("ITAG4.0_descriptions.txt", sep="\t", quote = "") # Load annotation data
#colnames(genes_annotation) <- c("ID", "Annotation")
#genes_annotation[,1] <- toupper(genes_annotation[,1]) # Change the ID to uppercase in the annotation data
#colnames(transcriptomics) <- gsub("\\.\\d+$", "", colnames(transcriptomics)) # Delete suffixes in the IDs
# Change the IDs for the gene names (column names):
#matched_annotations <- genes_annotation$Annotation[match(colnames(transcriptomics), genes_annotation$ID)]
#colnames(transcriptomics) <- ifelse(is.na(matched_annotations), colnames(transcriptomics), matched_annotations)
# Delete the parenthesis (example: "(AHRD V3.3 *** Q07459_SOLTU)")
#colnames(transcriptomics) <- gsub("\\s*\\(.*?\\)", "", colnames(transcriptomics))
# Save the resulting matrix as csv file
#write.csv(transcriptomics, file="transcriptomics.csv")

# Physiological data
pheno_archive <- "Phenotyping_data.xlsx" # Change with your archive name 
pheno_sheet <- "Sheet2"
pheno_data <- read_excel(pheno_archive, sheet = pheno_sheet)
pheno_data <- as.data.frame(pheno_data)
rownames(pheno_data) <- pheno_data[,1] # Set sample names as row names
pheno_data <- pheno_data[,-(1:3)]

damage_archive <- "Mite_damage.xlsx" # Change with your archive name 
damage_sheet <- "Sheet2"
Mite_damage <- read_excel(damage_archive, sheet = damage_sheet)
Mite_damage <- Mite_damage[,-(1:3)]

pheno_data$damage <- Mite_damage # 
pheno_data <- as.matrix(pheno_data)

# Define Y (factor that indicate the class of each sample)
anotacio_mostres <- read.csv("anotacio_mostres.csv", header=TRUE, sep=";")
Y <- as.factor(anotacio_mostres$Name)
Y <- gsub("_\\d+", "", Y)
table(Y)

# Delete variables with 0 variance
pheno_data <- pheno_data[, apply(pheno_data, 2, var) > 0, drop=FALSE]
transcriptomics <- transcriptomics[, apply(transcriptomics, 2, var) > 0, drop=FALSE]
metabolomics <- metabolomics[, apply(metabolomics, 2, var) > 0, drop=FALSE]

# Create a list with all the blocks 
X = list(physiological = pheno_data,
         transcriptomics = transcriptomics, 
         metabolomics = metabolomics)
lapply(X, dim)

# Assure that all the blocks are numeric
X$physiological <- apply(X$physiological, 2, as.numeric)
X$transcriptomics <- apply(X$transcriptomics, 2, as.numeric)
X$metabolomics <- apply(X$metabolomics, 2, as.numeric)


# ------------------------- Create the design matrix -------------------------------------
# To avoid "error: Unique indentifier is needed" in pls function
colnames(X$physiological) <- make.unique(colnames(X$physiological))
colnames(X$transcriptomics) <- make.unique(colnames(X$transcriptomics))
colnames(X$metabolomics) <- make.unique(colnames(X$metabolomics))

# Examine the correlation between blocks
pls1 <- pls(X$physiological, X$transcriptomics, ncomp = 1)
print("correlation of physiological and transcriptomics data")
cor(pls1$variates$X, pls1$variates$Y)

pls2 <- pls(X$physiological, X$metabolomics, ncomp = 1)
print("correlation of physiological and metabolomics data")
cor(pls2$variates$X, pls2$variates$Y)

pls3 <- pls(X$transcriptomics, X$metabolomics, ncomp = 1)
print("correlation of transcriptomics and metabolomics data")
cor(pls3$variates$X, pls3$variates$Y)

# Create the design matrix
# The choice of the design can be based on prior knowledge(the blocks you expect to be correlated) or data-driven
# based on the PLS analysis to examine the correlation between the blocks
weight <- 0.8 # Change accordingly with the results of the correlation or your interest
design = matrix(weight, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0 # Set the diagonal values to 0
design 

# ------------------------- Number of components ---------------------------------------
K <- length(table(Y)) # Number of classes in Y
# Create a basic model
basic.diablo.model = block.splsda(X, Y, 
                                  ncomp = 5, # You can change it to K-1
                                  design = design, 
                                  near.zero.var = TRUE) 

M <- 3 # Should be samples/M > 6
# Examine the performance of the model
perf.diablo = perf(basic.diablo.model, 
                   validation = 'Mfold', 
                   folds = M, 
                   nrepeat = 10)

# If you have a low number of samples (<10), change validation method to 'loo'
# perf.diablo = perf(basic.diablo.model, validation = "loo") 

# plot output of tuning
pdf(file="plot_perf_basic.diablo.model.pdf", width=10, height=10)
plot(perf.diablo) 
dev.off()
#Lists the different types of error rates
perf.diablo$error.rate

# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 

# Alternatively, you can manually select the number of components as K - 1
# ncomp <- K - 1

# ------------------------- Number of variables ---------------------------------------
# Set grid of values for each component to test (change accordingly to your data)
test.keepX = list(physiological = c(1, 2, 3),
                  transcriptomics = c(5:9, seq(100,1000, 50)),
                  metabolomics = seq(5, 30, 5))

# It may take some time to run, parallelisation is possible as follows:
#library(BiocParallel)
#BPPARAM <- MulticoreParam(workers = 4) # Change accorgingly to your CPU

# Run the feature selection tuning (it is advised to use nrepeat=10-50) (add BPPARAM for parallelisation)
tune = tune.block.splsda(X, Y, ncomp = ncomp, 
                         test.keepX = test.keepX, 
                         design = design,
                         validation = 'Mfold', folds= 3, 
                         nrepeat = 30,
                         progressBar = TRUE,
                         BPPARAM = SerialParam(),
                         dist = "centroids.dist")

# set the optimal values of features to retain
list.keepX = tune$choice.keepX 
list.keepX

# ------------------------- Final model ---------------------------------------
final.diablo.model <- block.splsda(X, Y, ncomp = ncomp,
                                   keepX = list.keepX, design = design)
# design matrix for the final model
final.diablo.model$design 

# Extract the selected variables for the first component (change 'comp' to other components)
physiological_selected <- selectVar(final.diablo.model, block = 'physiological', comp = 1)$physiological$name 
genes_selected <- selectVar(final.diablo.model, block = 'transcriptomics', comp = 1)$transcriptomics$name 
metabolites_selected <- selectVar(final.diablo.model, block = 'metabolomics', comp = 1)$metabolomics$name 

# ------------------------- Samples plots ---------------------------------------
# Save plots in pdf
pdf(file="DIABLO_sampleplots.pdf", width=6, height=6)

# check if the correlations between components were maximised as specified in the design matrix
plotDiablo(final.diablo.model, 
           ncomp = 1)

# understand the information extracted from each dataset and its discriminative ability.
plotIndiv(final.diablo.model, 
          ind.names = FALSE, 
          legend = TRUE, 
          title = 'DIABLO Sample Plots')

# highlight the agreement between all datasets at the sample level
plotArrow(final.diablo.model, 
          ind.names = FALSE, 
          legend = TRUE, 
          title = 'DIABLO Arrow plot')
dev.off()

# ------------------------- Variable plots ---------------------------------------
# Save plots in pdf
pdf(file="DIABLO_variableplots.pdf", width=9, height=9)

# highlights the contribution of each selected variable to each component
plotVar(final.diablo.model, 
        var.names = FALSE, 
        style = 'graphics', 
        legend = TRUE,
        pch = c(16, 17, 15), 
        cex = c(2,2,2), 
        col = c('darkolivegreen1', 'plum', 'cadetblue3')) # Change the colors of the blocks

# represents the correlations between variables of different types
circosPlot(final.diablo.model,
           comp = c(1:3), # select the components to represent
           cutoff = 0.9, # select the threshold
           line = TRUE,
           color.blocks= c('darkolivegreen1', 'plum', 'cadetblue3'),  # Change the colors of the blocks
           color.cor = c("chocolate3","grey20"), # Change the color of the correlations
           size.labels = 1.5)

# represent the multi-omics molecular signature expression for each sample
cimDiablo(final.diablo.model, 
          color.blocks = c('darkolivegreen1', 'plum', 'cadetblue3'), # Change the colors of the blocks
          comp = c(1:3), # select the components to represent
          margin=c(8,20), 
          legend.position = "right")

dev.off()


# Save the network to visualize in Cytoscape
my.network = network(final.diablo.model, 
                     blocks = c(1,2,3),
                     color.node = c('darkolivegreen1', 'plum', 'cadetblue3'), 
                     cutoff = 0.9, # Select the threshold
                     save = 'jpeg', 
                     name.save = 'diablo_network')
write_graph(my.network$gR, 
            file = "Network_09.gml", # Name of the archive generated
            format = "gml")
