
# Title: DIABLO step-by-step
# Created: November 2024
# Author: Paula Solé

# ---------------------------- Loading packages and data ------------------------------------

# Load packages
library(mixOmics)
library(gridExtra)
library(igraph)

# Clear your workplace
rm(list = ls())

# Change working directory
directory <- '/Volumes/ftp/Paula Sole/processed_data' # Change with your working directory
setwd(directory)

## ------- If your omics data is in the same dataset, follow this steps ----------

# Load data (samples in rows and features in columns)
data_archive <- "PA_dataframe.txt" # Change with your archive name
data <- read.table(data_archive, header=TRUE, sep="\t", row.names=1)
dim(data)

# Define Y (factor that indicate the class of each sample)
groups <- rownames(data)
Y <- as.factor(groups)
summary(Y)

# Delete columns with NA
columns_NA <- colSums(is.na(data)) > 0
index_columns_NA <- which(columns_NA)
print(index_columns_NA)
data <- data[, -index_columns_NA]

# Subset data to create the different blocks
physiological <- data[, (1:3)]
transcriptomics <- data[, grepl("SOLYC", colnames(data))]
metabolomics <- data[, grepl("met", colnames(data))]

## ------- If your omics data is in separates datasets, follow this steps ----------

# Metabolomics data
met_archive <- "Annotated_Met_QP.xlsx" # Change with your archive name 
met_sheet <- "Sheet1" 
metabolomics <-  read_excel(met_archive, sheet = met_sheet)
dim(metabolomics) # Check the dimension of the dat

# Transcriptomics data
trans_archive <- "processed_data.csv" # Change with your archive name 
transcriptomics <- read.csv(trans_archive, header=TRUE, sep=",")
dim(transcriptomics) # Check the dimension of the data

# Physiological data
physio_archive <- "Phenotyping_data.xlsx" # Change with your archive name 
physio_sheet <- "Sheet2"
physiological <- read_excel(physio_archive, sheet = physio_sheet)
dim(physiological) # Check the dimension of the data

# Define Y (factor that indicate the class of each sample)
groups_archive <- "anotacio_mostres.csv" # Change with your archive name
groups <- read.csv(groups_archive, header=TRUE, sep=";")
Y <- as.factor(anotacio_mostres$Name)
table(Y)

## ----------------------------- input data ----------------------------------- 

# Delete variables with 0 variance
physiological <- physiological[, apply(physiological, 2, var) > 0, drop=FALSE]
transcriptomics <- transcriptomics[, apply(transcriptomics, 2, var) > 0, drop=FALSE]
metabolomics <- metabolomics[, apply(metabolomics, 2, var) > 0, drop=FALSE]

# Create a list with the 3 blocks
X = list(physiological = physiological,
         transcriptomics = transcriptomics, 
         metabolomics = metabolomics)
lapply(X, dim)
# Alternatively, load each block separately and create the list

# ------------------------- Create the design matrix -------------------------------------

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

# Create the desing matrix
weight <- 0.7 # Change accorgingly with the results of the correlation or your interest
design = matrix(weight, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0
design 

# ------------------------- Number of components ---------------------------------------

K <- length(table(Y)) # Number of classes in Y
basic.diablo.model = block.splsda(X, Y, ncomp = K-1 , design = design) 

M <- 10 # Should be samples/M > 6
perf.diablo.tcga = perf(basic.diablo.model, validation = 'Mfold', folds = 10, nrepeat = 10)

# If you have a low number of samples (<10), change validation method to 'loo'
perf.diablo = perf(basic.diablo.model, validation = "loo") 

# plot output of tuning
plot(perf.diablo) 
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
# library(BiocParallel)
# BPPARAM <- MulticoreParam(workers = 4) # Change accorgingly to your CPU

# Run the feature selection tuning (it is advised to use nrepeat=10-50) (add BPPARAM for parallelisation)
tune = tune.block.splsda(X, Y, ncomp = ncomp, 
                         test.keepX = test.keepX, design = design,
                         validation = 'Mfold', folds= 3, nrepeat = 1, 
                         # BPPARAM = BPPARAM,
                         dist = "centroids.dist")

# set the optimal values of features to retain
list.keepX = tune$choice.keepX 
list.keepX


# ------------------------- Final model ---------------------------------------

final.diablo.model <- block.splsda(X, Y, ncomp = ncomp,
                                   keepX = list.keepX, design = design)
# design matrix for the final model
final.diablo.model$design 

# Extract the selected variables
physiological_selected <- selectVar(final.diablo.model, block = 'physiological', comp = 1)$physiological$name 
genes_selected <- selectVar(final.diablo.model, block = 'transcriptomics', comp = 1)$transcriptomics$name 
metabolites_selected <- selectVar(final.diablo.model, block = 'metabolomics', comp = 1)$metabolomics$name 

# ------------------------- Samples plots ---------------------------------------

# Save plots in pdf
pdf(file="~/Desktop/DIABLO_sampleplots.pdf", width=10, height=10)

# check if the correlations between components were maximised as specified in the design matrix
plotDiablo(final.diablo.model, ncomp = 1)

# understand the information extracted from each dataset and its discriminative ability.
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots')

# highlight the agreement between all datasets at the sample level
plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO')
dev.off()

# ------------------------- Variable plots ---------------------------------------

# Save plots in pdf
pdf(file="~/Desktop/DIABLO_variableplots.pdf", width=10, height=10)

# highlights the contribution of each selected variable to each component
plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'lightgreen', 'cadetblue2'))

# represents the correlations between variables of different types
circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'lightgreen', 'cadetblue2'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

# represent the multi-omics molecular signature expression for each sample
cimDiablo(final.diablo.model, 
          color.blocks = c('darkorchid', 'lightgreen', 'cadetblue2'), # Change the colors of the blocks
          comp = c(1:3), # select the components to represent
          margin=c(8,20), 
          legend.position = "right")

# visualise the correlations between the different types of variables
# X11() #Opens a new window
network(final.diablo.model, blocks = c(1,2,3),
        color.node = c('darkorchid', 'lightgreen', 'cadetblue2'), cutoff = 0.4,
        save = 'jpeg', name.save = 'diablo')

dev.off()

# Save the network to visualise in cytoscape
my.network = network(final.diablo.model, 
                     blocks = c(1,2,3),
                     color.node = c('darkorchid', 'lightgreen', 'cadetblue2'), 
                     cutoff = 0.4,
                     save = 'jpeg', 
                     name.save = 'diablo_network')
write_graph(my.network$gR, file = "Network_network.gml", format = "gml")
