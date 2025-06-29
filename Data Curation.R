#Loading the libraries required for downloading data
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)

#Check available data types
data_types <- curatedTCGAData(
    diseaseCode = "KIRC", assays = "*", version = "2.0.1")

##Dataset for Knowledge Guided Multi-Omic Network Inference (KiMONo)
#Downloading data as an MAE object
(kircmae <- curatedTCGAData(
       "KIRC", c("RNASeq2GeneNorm", "miRNASeqGene", "Methylation_methyl450"), version = "2.0.1", dry.run = FALSE
   ))

#Sample Types in the Dataset
sampleTables(kircmae)

#Split the MAE object into normal and tumor types
kirc_tumor <- TCGAsplitAssays(kircmae, "01")
kirc_normal <- TCGAsplitAssays(kircmae, "11")

#Subset the common samples
random_samples <- sample(names(which(table(sampleMap(kirc_tumor)$primary)==3)),150)
MAE_tumor <- subsetByColData(kirc_tumor, random_samples)
random_samples <- sample(names(which(table(sampleMap(kirc_normal)$primary)==3)),24)
MAE_normal <- subsetByColData(kirc_normal, random_samples)

#Exporting the dataset
td <- "C:/Documents/ccRCC"
tempd <- file.path(td, "Kimono")
if (!dir.exists(tempd))
  dir.create(tempd)

exportClass(MAE_tumor, dir = tempd, fmt = "csv", ext = ".csv")
exportClass(MAE_normal, dir = tempd, fmt = "csv", ext = ".csv")

##Dataset for IntOMICS
#Downloading data as an MAE object
(kircmae <- curatedTCGAData(
  "KIRC", c("RNASeq2GeneNorm", "Methylation_methyl450", "GISTIC_AllByGene"), version = "2.0.1", dry.run = FALSE
))

#Select only primary tumor
kirc_ma<- TCGAprimaryTumors(kircmae)

#Select matching samples across the three modalities
random_samples <- sample(names(which(table(sampleMap(kirc_ma)$primary)==3)),313)
kirc_ma<- subsetByColData(kirc_ma, random_samples)

rownames(kirc_ma[["KIRC_GISTIC_AllByGene-20160128"]]) <- rowData(kirc_ma[["KIRC_GISTIC_AllByGene-20160128"]])[["Gene.Symbol"]]