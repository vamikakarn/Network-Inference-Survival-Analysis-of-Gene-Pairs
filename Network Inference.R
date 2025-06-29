#Loading required libraries
library(kimono)
library(igraph)
library(data.table)
library(oem)
library(foreach)
library(doSNOW)
library(dplyr)  
library(tidyverse)
library(ggplot2)
library(cowplot) 
library(DT)

##Run the pipeline for both normal and tumor samples separately
#Load the data exported from the curatedTCGAData
GeneExp <- read.csv("GeneExp.csv", header = TRUE)
miRNAExp <- read.csv("miRNAExp.csv", header = TRUE)
Methylation <- read.csv("Methylation.csv", header = TRUE)

#Convert into data tables
GeneExp <- as.data.table(GeneExp)
miRNAExp <- as.data.table(miRNAExp)
Methylation <- as.data.table(Methylation)

#Create the input data as a list of data tables
input_data <- list (
  'gene' = GeneExp,
  'miRNA' = miRNA,
  'methylation' = Methylation
)

#Load the interaction files obtained from the respective databases 
gene_gene <- read.csv("gene_gene.csv", header  = TRUE)
gene_miRNA <- read.csv("gene_miRNA.csv", header =  TRUE)

#Create the prior network
prior_network <- create_prior_network(rbind(gene_proteome,gene_gene)) 

#Plotting the prior network
vertex <- do.call(rbind,strsplit(V(prior_network)$name,split = '___'))

prior_network %>% plot(edge.curved=0,
                       main = 'Prior Network',
                       vertex.color = c("steel blue", "orange")[vertex[,1] %>% as.factor %>% as.numeric],
                       vertex.frame.color="white",
                       vertex.label = vertex[,2], 
                       vertex.label.color='black',
                       vertex.label.cex=.7,
                       layout=layout_randomly, rescale=F) 
legend(x=-1.5, y=-1.1, c("Genes","miRNA"), pch=21,
       col="#777777", pt.bg=c("steel blue", "orange"), pt.cex=2, cex=.8, bty="n", ncol=1)

#Call KiMONo for network inference using SGL
network <- kimono(input_data, prior_network ,core = 2, infer_missing_prior = TRUE)

#Export the network file with contains value (effect size), r-squared (model performance), and mse (error) for further analysis
write.csv(network, file = "network.csv")


                        
                       
                       