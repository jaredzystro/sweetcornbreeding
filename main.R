
### Workflow for data analysis and prediction of 
### Jared Zystro's se sweet corn breeding project

# Style guide:
# ALLCAPS - global variables
# UpperCamel - functions
# lower_underscore - variables
# Space between operators
# No space after function call

# PART A - for each trait

# 1. Initialize
#     Clean up
#     Set global variables
InitializeProject <- function() {
  
  rm(list = ls()) 
  
  WORKING_DIRECTORY <<- "~/Dropbox/PhD/SweetCornBreeding"
  PHENO_DATA_FILE <<- "./data/Data.csv"
  MARKER_DATA_FILE <<- "./data/data_parimp_37genos_dropminorhets_dropNAs_thinto250KB.hmp.txt"
  FUNCTIONS_SUBDIRECTORY <<- "/functions"
  
  setwd(WORKING_DIRECTORY)
  
}
InitializeProject()

# 2. Import / create data tables
#     Load in 2015 and 2016 phenotype data
LoadPhenoData <- function() {
  
  return(read.csv(PHENO_DATA_FILE))
  
}
pheno_data <- LoadPhenoData()

#     Load in inbred marker data
LoadMarkerData <- function() {
  
  return(read.delim(MARKER_DATA_FILE))
  
}
marker_data <- LoadMarkerData()

#     Create pedigree table

# 3. Format data
#     Set factors in phenotype data
AddPhenoDataFactors <- function (pheno_data) {
  
pheno_data$Geno <- factor(pheno_data$Geno)
pheno_data$LocBlock <- paste0(pheno_data$Location,pheno_data$Block)
pheno_data$Block <- factor(pheno_data$Block)
pheno_data$LocRep <- paste0(pheno_data$Location,pheno_data$Rep)
pheno_data$LocRep <- factor(pheno_data$LocRep)
pheno_data$Set <- factor(pheno_data$DIISet)
pheno_data$MaleGeno <- factor(pheno_data$MaleGeno)
pheno_data$FemaleGeno <- factor(pheno_data$FemaleGeno)
pheno_data$Geno <- factor(pheno_data$Geno)

return (pheno_data)

}
pheno_data <- AddPhenoDataFactors(pheno_data)

#     Format raw marker data
FormatMarkerData <- function (marker_data) {
  
  # Convert from HapMap format into 0,1,2 format
  # Code from the NAM package vignette: 
  # ftp://ftp.math.ethz.ch/sfs/CRAN/web/packages/NAM/vignettes/vignette.html
  G <- marker_data
  AA = as.character(G[,2]); 
  AA = gsub('/','',AA); 
  Gen = t(G[,-c(1:11)]); 
  n=nrow(Gen); 
  m=ncol(Gen); 
  gen = matrix(NA,n,m); 
  for(i in 1:m){ 
    A1 = strsplit(AA[i],'')[[1]][1]; 
    A2 = strsplit(AA[i],'')[[1]][2];
    BB = paste(A1,A1,sep=''); 
    Bb = paste(A1,A2,sep=''); 
    bb = paste(A2,A2,sep=''); 
    M = as.character(Gen[,i]); 
    gen[M==BB,i]=2; 
    gen[M==Bb,i]=1;
    gen[M==bb,i]=0}; 
  colnames(gen)=paste(G[,3],G[,4],G[,2],sep='.')
  gen2<-gen
  
  # Clean up names
  rownames(gen2) <- gsub("(.*?)\\.(.*)","\\1",rownames(Gen))
  rownames(gen2)[rownames(gen2)=="We1095047"]<-"We1095407-4"
  rownames(gen2)[rownames(gen2)=="We1095407"]<-"We1095407-2"
  rownames(gen2)[rownames(gen2)=="We09423"]<-"We09423a"
  rownames(gen2)[rownames(gen2)=="We07812"]<-"Wu07812"
  rownames(gen2)[rownames(gen2)=="We10802R"]<-"Wu10802R"
  rownames(gen2)[rownames(gen2)=="We10803R"]<-"Wu10803R"
  
  marker_data <- gen2
  
  return (marker_data)
  
}
marker_data <- FormatMarkerData(marker_data)

#     Create marker data for hybrids from inbred marker data

# 4. Quality control
#     Check for outliers
#     Test for normality
#     Test for equal variance

# 5. Test for experient-wide GxE
#     ANOVA
#     Spearman rank correlation

# 6. Test to see if means should be adjusted based on checks
#     ANOVAs of checks

# 7. ANOVAs
#     Inbred ANOVAs
#     Hybrid ANOVAs
 
# 8. Collect phenotype values
#     Inbred means
#     Tested hybrid means
#     Inbred GCAs
#     Tested hybrid SCAs

# 9. Predict untested hybrid means
#     Classical (u + GCAf + GCAm)
#     GBLUP
#     Cross validiate GBLUP
#     Classical vs. GBLUP comparison

# 10. Predict synthetics
#     Classical (two step)
#     GBLUP (two step)
#     GLBUP (one step)

# PART B - for all traits

# 11. Merge trait data

# 12. Filter and sort hybrids and synthetics based on index

# 13. Create tables and charts
#      QC info - QQplots, boxplots, etc
#      ANOVA tables
#      Spearkman correlations
#      Inbred means
#      Tested hybrid means
#      Inbred GCAs
#      Tested hybrid SCAs
#      Predicted hybrid means (Classical and GBLUP)
#      GBLUP cross-validation 
#      Predicted synthetic means (filtered)






