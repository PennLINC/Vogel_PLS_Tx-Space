# Script to process GTEx data (v7)
# Gokul Ramaswami 5-10-2018, adapted for use on GTEx v8 by Jakob Seidlitz

# set up R
rm(list=ls())
options(stringsAsFactors=FALSE)

# Load libraries
library(WGCNA)
library(SNFtool)
library(flashClust)
library(CePa)
library(plyr)
#enableWGCNAThreads(nThreads=8); # when requesting 8 CPUs

######################### FUNCTIONS ##########################

## Covariate and expression PC relationship functions

# Run PCA funtion
Run_PCA <- function (data, nPC = 5) {
  pCdat <- prcomp(t(data), center=FALSE,scale=FALSE);
  topPCs <- pCdat$x[,1:nPC];
  # Calculate variance explained by each PC
  varExp <- (pCdat$sdev)^2 / sum(pCdat$sdev^2)
  topVar <- varExp[1:nPC]
  colnames(topPCs) <- paste("Data\n", colnames(topPCs), " (", signif(100 * topVar[1:nPC], 2), "%)", sep = "")
  return(topPCs)
}

# Correlation panel function
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if ((class(x) == "numeric" | class(x) == "integer") & (class(y) == "numeric" | class(y) == "integer")) {
    r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1.2/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# Plot correlation panel function
Plot_Corr_Matrix <- function (pairsDat, topPCs, title, outPath) {
  pdf(outPath, height = 20, width = 24)
  pairs(cbind(pairsDat, topPCs), pch = 19, upper.panel = panel.cor, main = paste("Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values","\n", title, sep = ""),cex.labels=1)
  dev.off()
}

# Function to perform linear regression covariate correction of Expression data
regressCovariates_mRNA_linearModel <- function(datExpr,datMeta) {
  
  # choose covariates to regress
  RIN <- as.numeric(datMeta[,"RIN"])
  age <- as.numeric(datMeta[,"AGE"])
  sex <- as.numeric(as.factor(datMeta[,"SEX"]))-1
  seqpc1 <- as.numeric(datMeta[,"SeqPC1"])
  seqpc2 <- as.numeric(datMeta[,"SeqPC2"])
  seqpc3 <- as.numeric(datMeta[,"SeqPC3"])
  seqpc4 <- as.numeric(datMeta[,"SeqPC4"])
  seqpc5 <- as.numeric(datMeta[,"SeqPC5"])
  # for brain region - separate into 2 surrogate variables: BA24 and BA9
  BA24 <- as.numeric(datMeta[,"SMTSD"]=="Brain - Anterior cingulate cortex (BA24)")
  BA9 <- as.numeric(datMeta[,"SMTSD"]=="Brain - Frontal Cortex (BA9)")
  
  regvars <- as.data.frame(cbind(RIN,age,sex,seqpc1,seqpc2,seqpc3,seqpc4,seqpc5))
  
  ## Run the regression and make the adjusted matrix via matrix multiplication                                                                                                                                       
  
  X <- model.matrix(as.formula(paste0("~",paste(colnames(regvars),collapse="+"))), data = regvars)
  rownames(X) <- rownames(datMeta)
  Y <- datExpr
  
  beta <- (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  
  # regress out beta values from data
  to_regress <- (as.matrix(X[,2:(ncol(regvars)+1)]) %*% (as.matrix(beta[2:(ncol(regvars)+1),])))
  datExpr.reg <- datExpr - t(to_regress)
  
  ## This datExpr.reg is now a technical variable corrected matrix.
  rownames(datExpr.reg) <- rownames(datExpr)
  colnames(datExpr.reg) <- rownames(datMeta)
  
  datExpr.reg
  
}


################################################################

# define directories
dataDir <- "~/Desktop/Education_GTEX_WGCNA/data/"
figureDir <- "~/Desktop/Education_GTEX_WGCNA/figures/"
if (!dir.exists(figureDir)) {
  dir.create(figureDir)
}
processedDataDir <- "~/Desktop/Education_GTEX_WGCNA/data/"
if (!dir.exists(processedDataDir)) {
  dir.create(processedDataDir)
}


###### Load in GTEX Data and pre-processing ######


# Load in GTEX data
if (!file.exists(paste0(processedDataDir,"InputData.RData"))) { 
  gene_tpm <- read.gct(paste0(dataDir,'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'))
  meta_attr <- read.table(paste0(dataDir,'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'), sep = "\t",header=TRUE,quote="") # this was very annoying to figure out why the file kept getting truncated - need the quote="" to avoid premature EOF
  #gene_tpm <- read.gct(paste0(dataDir,'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct'))
  #meta_attr <- read.table(paste0(dataDir,'GTEx_v7_Annotations_SampleAttributesDS.txt'), sep = "\t",header=TRUE,quote="") # this was very annoying to figure out why the file kept getting truncated - need the quote="" to avoid premature EOF
  # get cortex sample only
  #meta_attr_ctx <- meta_attr[meta_attr[,"SMTSD"] %in% c("Brain - Cortex","Brain - Frontal Cortex (BA9)","Brain - Anterior cingulate cortex (BA24)"),]
  meta_attr_ctx <- meta_attr[meta_attr[,"SMTSD"] %in% c(levels(as.factor(meta_attr$SMTSD))[grep("Brain",levels(as.factor(meta_attr$SMTSD)))[-12]]),]
  # change annoying gtex naming conventions so sample attributes match)
  colnames(gene_tpm) <- gsub("\\.","-",colnames(gene_tpm))
  gene_tpm_ctx <- gene_tpm[,colnames(gene_tpm) %in% meta_attr_ctx$SAMPID]
  rownames(meta_attr_ctx) <- meta_attr_ctx$SAMPID
  meta_attr_ctx = meta_attr_ctx[colnames(gene_tpm_ctx),]
  save(gene_tpm_ctx,meta_attr_ctx,file=paste0(processedDataDir,"InputData.RData"))
} else {
  load(paste0(processedDataDir,"InputData.RData"))
}

# remove zero variance genes
gene_tpm_ctx2 <- gene_tpm_ctx[apply(gene_tpm_ctx,1,var)>0,]
# log2 transform expression to stabilize variance
gene_tpm_ctx3 <- log2(gene_tpm_ctx2+0.001)
# check median expression of these genes
expr_med <- apply(gene_tpm_ctx3,1,median)
# filter to those genes with median log2(expression) > 0 --- remove lowly expressed genes in cortex
gene_tpm_ctx4 <- gene_tpm_ctx3[expr_med > 0,]
# fix gene names - remove trailing .number
rownames(gene_tpm_ctx4) <- gsub("\\.\\d+$","",rownames(gene_tpm_ctx4),perl=TRUE)

# Get sequencing statsitics to calculate sequencing PCS
# columns 19-63
seq_param <- meta_attr_ctx[,19:63] 
# remove NA columns
seq_param <- seq_param[,colSums(is.na(seq_param))<nrow(seq_param)]
# remove zero variance columns
seq_param <- seq_param[,apply(seq_param,2,sd) > 0]

# log transform the LARGE number columns: 5, 9, 11-15, 19, 21-23, 25-26, 28)
seq_param_log <- seq_param
seq_param_log[,c(5,9,11:15,19,21:23,25:26,28)] = log10(seq_param_log[,c(5,9,11:15,19,21:23,25,26,28)])
seq_param_log_scaled <- scale(seq_param_log)
seq.pca <- prcomp(seq_param_log_scaled)

# get top 5 PC's
num_topPCs <- 5 # Number of top seqStatPCs to test
topPC.datSeq <- seq.pca$x[,1:num_topPCs];
colnames(topPC.datSeq) <- paste0("SeqPC",c(1:num_topPCs)) ## Compute sequencing stats PCA
rownames(topPC.datSeq) <- meta_attr_ctx$SAMPID

## get subject info
subject_meta <- read.table(paste0(dataDir,'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'), sep ='\t', header = TRUE)
rownames(subject_meta) <- subject_meta$SUBJID
# in the absence of the non-public data take the mean age for each bin
subject_meta$AGE <- as.numeric(as.character(revalue(subject_meta$AGE, c("20-29"="25", "30-39"="35",
                                                                        "40-49"="45", "50-59"="55",
                                                                        "60-69"="65", "70-79"="75"))))
# add subjectID to metadata
ID_SPLIT <- vector(length=nrow(meta_attr_ctx))
for ( i in 1:length(ID_SPLIT)) {
  splitName <- unlist(strsplit(rownames(meta_attr_ctx)[i],"-"))
  ID_SPLIT[i] = paste(splitName[1],splitName[2],sep="-")
}
meta_attr_ctx$SUBJID <- ID_SPLIT
# Add Age and Sex to metadata
meta_attr_ctx$AGE <- subject_meta[meta_attr_ctx$SUBJID,"AGE"]
meta_attr_ctx$SEX <- subject_meta[meta_attr_ctx$SUBJID,"SEX"]

# Subset to covariates of interest
ctx_covariates <- cbind(meta_attr_ctx[,c("SUBJID","SMRIN","SMNABTCH","SMGEBTCH","AGE","SEX","SMTSD")],topPC.datSeq)
colnames(ctx_covariates)[c(2:4)] <- c("RIN","BATCH","GBATCH")
ctx_covariates$SEX = as.factor(ctx_covariates$SEX)
ctx_covariates$AGE = as.numeric(ctx_covariates$AGE)
ctx_covariates$BATCH = as.factor(ctx_covariates$BATCH)
ctx_covariates$GBATCH = as.factor(ctx_covariates$GBATCH)
ctx_covariates$SMTSD = as.factor(ctx_covariates$SMTSD)

###### Regression of covariates and outlier removal ######


## Check relationship between top gene expression PCs and covariates
# Standard normalize unregressed data
gene_tpm_ctx4.stdNorm <- t(standardNormalization(t(gene_tpm_ctx4)))
write.csv(gene_tpm_ctx4.stdNorm,"~/Desktop/Education_GTEX_WGCNA/processed_data/GTEx_brain_norm_noregress.csv")
# Set up covariate matrix
# pairsDat_mRNA <- ctx_covariates[,c("RIN","AGE","SEX","SMTSD","SeqPC1","SeqPC2","SeqPC3","SeqPC4","SeqPC5")]
# # Unnormalized Expression Data
# topPCs_mRNA_unNormalized <- Run_PCA(gene_tpm_ctx4.stdNorm)
# # Plot_Corr_Matrix(pairsDat_mRNA, topPCs_mRNA_unNormalized, title = "mRNA unregressed", outPath = paste0(figureDir, "Expression_Covariates_Unregressed.pdf"))
# 
# # Regress out technical covariates for expression data
# gene_tpm_ctx4.reg <- regressCovariates_mRNA_linearModel(gene_tpm_ctx4,ctx_covariates)
# # Standard normalize regressed data
# gene_tpm_ctx4.reg.stdNorm <- t(standardNormalization(t(gene_tpm_ctx4.reg)))
# # Unnormalized Expression Data
# topPCs_mRNA_regressed <- Run_PCA(gene_tpm_ctx4.reg.stdNorm)
# # Plot_Corr_Matrix(pairsDat_mRNA, topPCs_mRNA_regressed, title = "mRNA Covar Regressed", outPath = paste0(figureDir, "Expression_Covariates_CovarRegressed.pdf"))
# 
# # Remove Outlier samples (check using sample-sample connectivity)
# sdout <- 2; normadj <- (0.5+0.5*bicor(gene_tpm_ctx4.reg, use='pairwise.complete.obs'))^2
# netsummary <- fundamentalNetworkConcepts(normadj); 
# K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
# C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
# outliers_CTX <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
# print(paste("There are ",sum(outliers_CTX)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(gene_tpm_ctx4.reg)[outliers_CTX]); print(table(outliers_CTX))
# # pdf(paste0(figureDir,"Outliers_ZConnect.pdf"))
# # plot(Z.K, col = "black", pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
# # abline(h=sdout*-1, lty=2)
# # dev.off()
# 
# # remove outliers
# datExpr = gene_tpm_ctx4.reg[,!outliers_CTX]
# datMeta = ctx_covariates[!outliers_CTX,]
