#For the 3 target genes and their respective predictors IntOMICS was applied to run MCMC Simulation

# Load required libraries
library(IntOMICS)
suppressPackageStartupMessages({
  library(curatedTCGAData)
  library(TCGAutils)
  library(bnlearn)
  library(bnstruct)
  library(matrixStats)
  library(parallel)
  library(RColorBrewer)
  library(bestNormalize)
  library(igraph)
  library(gplots)
  library(methods)
  library(ggraph)
  library(ggplot2)})

#Load the annotation and gene annotation file. Annotation is a named list of corresponding probe IDs of associated genes. Gene Annotation is a data frame consisting of entrez IDs of the corresponding genes.
annot <- list(GeneName = c("Enter Probe IDs"))
gene_annot <- data.frame(entrezID = c("Enter entrez IDS"), gene_symbol = c("Enter gene symbols") )

data(list=c("gene_annot","annot"))
rowselect <- list(gene_annot$gene_symbol, gene_annot$gene_symbol, 
                  unlist(annot))

#Select only genes of interest
names(rowselect) <- names(kirc_ma)
omics <- kirc_ma[rowselect, , ]
names(omics) <- c("cnv","ge","meth")

#Load the prior network of gene-gene interactions obtained from respective database containing the columns src_entrez, dest_entrez, and edge_type
PK <- read.csv("gene_gene.csv", header  = TRUE)

#This framework also requires layers_def and TFtarg_mat. TFtarg_mat is a logical matrix containing known transcription factors and their targets.layers_def is a data frame. Both are already present in the package, however, can use others as required according to the data.
data("layers_def")
data("TFtarg_mat")

#Data Pre-processing to compute biological prior matrix
OMICS_mod_res <- omics_module(omics = omics, 
                              PK = PK, 
                              layers_def = layers_def, 
                              TFtargs = TFtarg_mat,
                              annot = annot, 
                              gene_annot = gene_annot,
                              lm_METH = TRUE,
                              r_squared_thres = 0.3,
                              p_val_thres = 0.1)

#MCMC Simulation to compute empirical biological knowledge matrix
if(interactive())
{
  BN_mod_res_sparse <- bn_module(burn_in = 100.000, 
                                 thin = 500, 
                                 OMICS_mod_res = OMICS_mod_res,
                                 minseglen = 50.000)
  
#Generate MCMC diagnostics graphs
trace_plots(mcmc_res = BN_mod_res,
              burn_in = 10000,
              thin = 500, 
              edge_freq_thres = 0.5)

#Generate network structures and heatmap
res_weighted <- edge_weights(mcmc_res = BN_mod_res, 
                             burn_in = 10000, 
                             thin = 500, 
                             edge_freq_thres = 0.5)

weighted_net_res <- weighted_net(cpdag_weights = res_weighted,
                                 gene_annot = gene_annot, 
                                 PK = PK, 
                                 OMICS_mod_res = OMICS_mod_res, 
                                 gene_ID = "gene_symbol", 
                                 TFtargs = TFtarg_mat, 
                                 B_prior_mat_weighted = B_prior_mat_weighted(BN_mod_res))
ggraph_weighted_net(net = weighted_net_res, 
                    node_size = 10, 
                    node_label_size = 4, 
                    edge_label_size = 4)

weighted_net_res <- weighted_net(cpdag_weights = res_weighted,
                                 gene_annot = gene_annot, 
                                 PK = PK, 
                                 OMICS_mod_res = OMICS_mod_res, 
                                 gene_ID = "gene_symbol", 
                                 edge_weights = "empB",
                                 TFtargs = TFtarg_mat, 
                                 B_prior_mat_weighted = B_prior_mat_weighted(BN_mod_res))#Graphs with edge weights as empirical prior knowledge inferred by IntOMICS.

ggraph_weighted_net(net = weighted_net_res)

#Heatmap for the difference between empirical biological knowledge and biological prior knowledge in gene-gene interactions
emp_b_heatmap(mcmc_res = BN_mod_res, 
              OMICS_mod_res = OMICS_mod_res, 
              gene_annot = gene_annot, 
              TFtargs = TFtarg_mat)


#Density of the edge weights inferred by IntOMICS
res_weighted <- edge_weights(mcmc_res = BN_mod_res, 
                             burn_in = 10000, 
                             thin = 500,
                             edge_freq_thres = NULL)

weighted_net_res <- weighted_net(cpdag_weights = res_weighted,
                                 gene_annot = gene_annot, 
                                 PK = PK, 
                                 OMICS_mod_res = OMICS_mod_res, 
                                 gene_ID = "gene_symbol", 
                                 TFtargs = TFtarg_mat, 
                                 B_prior_mat_weighted = B_prior_mat_weighted(BN_mod_res))

dens_edge_weights(weighted_net_res)


