
### Workflow for data analysis and prediction of 
### Jared Zystro's se sweet corn breeding project

# Style guide:
# ALLCAPS - global variables
# UpperCamel - functions
# lower_underscore - variables

# PART A - for each trait

# 1. Initialize
#     Clean up
#     Set global variables
InitializeProject <- function() {
  
  rm (list=ls()) 
  
  WORKING_DIRECTORY <<- "~/Dropbox/PhD/SweetCornBreeding"
  PHENO_DATA_FILE <<- "/data/Data.csv"
  MARKER_DATA_FILE <<- "/data/data_parimp_37genos_dropminorhets_dropNAs_thinto250KB.hmp.txt"
  FUNCTIONS_SUBDIRECTORY <<- "/functions"
  
  setwd(WORKING_DIRECTORY)
  
}
InitializeProject

# 2. Import / create data tables
#     Load in 2015 and 2016 phenotype data
LoadPhenoData <- function() {
  
  return(read.csv(PHENO_DATA_FILE))
  
}
pheno_data<-LoadPhenoData

#     Load in inbred marker data
LoadMarkerData <- function() {
  
  
  
}
marker_data <- LoadMarkerData

#     Create pedigree table

# 3. Format data
#     Set factors in phenotype data
#     Format raw marker data
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






